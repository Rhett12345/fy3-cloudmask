"""Command-line interface for cloud mask processing.

Usage:
    # Process single orbit
    python -m fy3_cloudmask process --config config/default.yaml --l1b data/L1b.HDF --geo data/GEO.HDF --output output/

    # Process batch
    python -m fy3_cloudmask batch --config config/default.yaml --start 2023-01-01 --end 2023-01-31 --data data/ --output output/

    # Convert thresholds
    python -m fy3_cloudmask convert-thresholds --input coeff/fylat_thresholds.mersi.ii3d.v8 --output config/thresholds/mersi_ii3d_v8.yaml
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

try:
    import click
    HAS_CLICK = True
except ImportError:
    HAS_CLICK = False

logger = logging.getLogger(__name__)


def setup_logging(verbose: bool = False) -> None:
    """Setup logging configuration.

    Args:
        verbose: Enable verbose logging.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )


if HAS_CLICK:
    @click.group()
    @click.option('--verbose', '-v', is_flag=True, help='Enable verbose logging')
    @click.pass_context
    def cli(ctx, verbose):
        """FY-3D MERSI-II Cloud Mask Processing Tool."""
        setup_logging(verbose)
        ctx.ensure_object(dict)
        ctx.obj['verbose'] = verbose

    @cli.command()
    @click.option('--config', '-c', required=True, help='Path to YAML config file')
    @click.option('--l1b', required=True, help='L1b HDF5 file path')
    @click.option('--geo', required=True, help='GEO HDF5 file path')
    @click.option('--output', '-o', required=True, help='Output directory')
    @click.option('--nwp1', help='First NWP GRIB file (optional)')
    @click.option('--nwp2', help='Second NWP GRIB file (optional)')
    @click.option('--oisst', help='OISST file (optional)')
    def process(config, l1b, geo, output, nwp1, nwp2, oisst):
        """Process a single FY-3D MERSI-II orbit."""
        from .pipeline import CloudMaskPipeline

        logger.info(f"Processing orbit: {l1b}")

        pipeline = CloudMaskPipeline(config)
        result = pipeline.process_orbit(
            l1b_path=l1b,
            geo_path=geo,
            output_dir=output,
            nwp_path1=nwp1,
            nwp_path2=nwp2,
            oisst_path=oisst,
        )

        if result.success:
            logger.info(f"Processing complete: {result.output_path}")
            logger.info(f"  Time: {result.processing_time:.1f}s")
            logger.info(f"  Cloudy: {result.n_cloudy} ({100*result.n_cloudy/result.n_pixels_processed:.1f}%)")
            sys.exit(0)
        else:
            logger.error(f"Processing failed: {result.error_message}")
            sys.exit(1)

    @cli.command()
    @click.option('--config', '-c', required=True, help='Path to YAML config file')
    @click.option('--start', required=True, help='Start date (YYYY-MM-DD)')
    @click.option('--end', required=True, help='End date (YYYY-MM-DD)')
    @click.option('--data', '-d', required=True, help='Data root directory')
    @click.option('--output', '-o', required=True, help='Output root directory')
    @click.option('--workers', '-n', default=1, help='Number of parallel workers')
    def batch(config, start, end, data, output, workers):
        """Process multiple orbits in batch mode."""
        from .pipeline import CloudMaskPipeline

        logger.info(f"Batch processing: {start} to {end}")

        pipeline = CloudMaskPipeline(config)
        results = pipeline.process_batch(
            start_date=start,
            end_date=end,
            data_root=data,
            output_root=output,
            n_workers=workers,
        )

        n_success = sum(1 for r in results if r.success)
        n_failed = sum(1 for r in results if not r.success)

        logger.info(f"Batch complete: {n_success} succeeded, {n_failed} failed")

        if n_failed > 0:
            for r in results:
                if not r.success:
                    logger.error(f"  Failed: {r.l1b_path} - {r.error_message}")

        sys.exit(0 if n_failed == 0 else 1)

    @cli.command('convert-thresholds')
    @click.option('--input', '-i', 'input_path', required=True, help='Input threshold file')
    @click.option('--output', '-o', 'output_path', required=True, help='Output YAML file')
    def convert_thresholds(input_path, output_path):
        """Convert Fortran threshold file to YAML format."""
        from .scripts.convert_thresholds import convert_threshold_file

        logger.info(f"Converting thresholds: {input_path} -> {output_path}")
        convert_threshold_file(input_path, output_path)
        logger.info("Conversion complete")

    def main():
        """Main entry point for CLI."""
        cli()

else:
    # Fallback if click is not installed
    def main():
        """Main entry point for CLI (click not available)."""
        print("Error: click package is required for CLI usage")
        print("Install with: pip install click")
        sys.exit(1)


if __name__ == '__main__':
    main()
