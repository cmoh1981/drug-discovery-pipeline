"""CLI entry point: python -m drugdiscovery"""

import sys

from drugdiscovery.config import parse_args
from drugdiscovery.pipeline import run_pipeline


def main() -> None:
    cfg = parse_args()
    try:
        run_pipeline(cfg)
    except KeyboardInterrupt:
        print("\nPipeline interrupted by user.")
        sys.exit(1)
    except Exception as exc:
        print(f"\nPipeline failed: {exc}")
        sys.exit(2)


if __name__ == "__main__":
    main()
