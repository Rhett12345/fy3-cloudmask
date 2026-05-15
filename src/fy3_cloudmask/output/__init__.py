"""Output module for writing cloud mask products."""

from .writer import write_cloud_mask, write_cloud_amount, write_combined_product
from .cloud_amount import compute_cloud_amount, compute_cloud_amount_with_coords

__all__ = [
    'write_cloud_mask',
    'write_cloud_amount',
    'write_combined_product',
    'compute_cloud_amount',
    'compute_cloud_amount_with_coords',
]
