"""Test fixtures for the service layer."""

from __future__ import annotations

import os

import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine, event
from sqlalchemy.orm import sessionmaker

# Override settings BEFORE importing service modules
os.environ["DATABASE_URL"] = "sqlite:///./test_service.db"
os.environ["SECRET_KEY"] = "test-secret-key-not-for-production"
os.environ["PIPELINE_OUTPUT_DIR"] = "test_job_results"

from service.app import create_app
from service.database import Base, get_db
from service.auth import hash_password
from service.models import User


# In-memory SQLite for tests
TEST_DATABASE_URL = "sqlite:///./test_service.db"
test_engine = create_engine(
    TEST_DATABASE_URL,
    connect_args={"check_same_thread": False},
)


@event.listens_for(test_engine, "connect")
def _set_sqlite_wal(dbapi_conn, connection_record):
    cursor = dbapi_conn.cursor()
    cursor.execute("PRAGMA journal_mode=WAL")
    cursor.execute("PRAGMA foreign_keys=ON")
    cursor.close()


TestSessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=test_engine)


@pytest.fixture(autouse=True)
def setup_test_db():
    """Create tables before each test, drop after."""
    Base.metadata.create_all(bind=test_engine)
    yield
    Base.metadata.drop_all(bind=test_engine)
    # Clean up test db file
    import pathlib
    db_file = pathlib.Path("test_service.db")
    if db_file.exists():
        try:
            db_file.unlink()
        except PermissionError:
            pass


def override_get_db():
    db = TestSessionLocal()
    try:
        yield db
    finally:
        db.close()


@pytest.fixture
def client():
    """FastAPI test client with overridden DB."""
    app = create_app()
    app.dependency_overrides[get_db] = override_get_db
    with TestClient(app) as c:
        yield c


@pytest.fixture
def test_db():
    """Direct DB session for test setup."""
    db = TestSessionLocal()
    try:
        yield db
    finally:
        db.close()


@pytest.fixture
def test_user(test_db):
    """Create a test user and return (user, plain_password)."""
    password = "testpassword123"
    user = User(
        email="test@example.com",
        hashed_password=hash_password(password),
        full_name="Test User",
        organization="Test Org",
    )
    test_db.add(user)
    test_db.commit()
    test_db.refresh(user)
    return user, password


@pytest.fixture
def auth_headers(client, test_user):
    """Get authorization headers for the test user."""
    user, password = test_user
    resp = client.post("/api/auth/login", json={
        "email": user.email,
        "password": password,
    })
    token = resp.json()["access_token"]
    return {"Authorization": f"Bearer {token}"}
