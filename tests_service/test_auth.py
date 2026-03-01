"""Tests for authentication endpoints."""

from __future__ import annotations


def test_register_success(client):
    resp = client.post("/api/auth/register", json={
        "email": "newuser@example.com",
        "password": "securepassword123",
        "full_name": "New User",
        "organization": "Acme Labs",
    })
    assert resp.status_code == 201
    data = resp.json()
    assert data["email"] == "newuser@example.com"
    assert data["full_name"] == "New User"
    assert "id" in data


def test_register_duplicate_email(client, test_user):
    user, _ = test_user
    resp = client.post("/api/auth/register", json={
        "email": user.email,
        "password": "anotherpassword123",
    })
    assert resp.status_code == 409


def test_register_short_password(client):
    resp = client.post("/api/auth/register", json={
        "email": "short@example.com",
        "password": "short",
    })
    assert resp.status_code == 422


def test_login_success(client, test_user):
    user, password = test_user
    resp = client.post("/api/auth/login", json={
        "email": user.email,
        "password": password,
    })
    assert resp.status_code == 200
    data = resp.json()
    assert "access_token" in data
    assert data["token_type"] == "bearer"


def test_login_wrong_password(client, test_user):
    user, _ = test_user
    resp = client.post("/api/auth/login", json={
        "email": user.email,
        "password": "wrongpassword123",
    })
    assert resp.status_code == 401


def test_login_nonexistent_user(client):
    resp = client.post("/api/auth/login", json={
        "email": "nobody@example.com",
        "password": "somepassword123",
    })
    assert resp.status_code == 401


def test_get_profile(client, auth_headers):
    resp = client.get("/api/auth/me", headers=auth_headers)
    assert resp.status_code == 200
    data = resp.json()
    assert data["email"] == "test@example.com"


def test_get_profile_unauthorized(client):
    resp = client.get("/api/auth/me")
    assert resp.status_code in (401, 403)


def test_get_profile_invalid_token(client):
    resp = client.get("/api/auth/me", headers={"Authorization": "Bearer invalid-token"})
    assert resp.status_code == 401
