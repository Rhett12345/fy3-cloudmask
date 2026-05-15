"""Allow running the package as a module.

Usage:
    python -m fy3_cloudmask process --config config/default.yaml --l1b data/L1b.HDF --geo data/GEO.HDF --output output/
"""

from .cli import main

if __name__ == '__main__':
    main()
