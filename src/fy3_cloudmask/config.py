"""Configuration management for the FY-3D cloud mask system.

Replaces the Fortran namelist (.nml) format with YAML configuration.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import yaml


@dataclass
class SensorConfig:
    """Sensor identification and geometry."""
    sensor_id: int = 21
    n_elem: int = 2048
    n_line: int = 2000


@dataclass
class PathConfig:
    """Input/output data paths."""
    code_root: str = ""
    l1b_data_dir: str = ""
    geo_data_dir: str = ""
    nwp_data_dir: str = ""
    oisst_data_dir: str = ""
    output_dir: str = ""
    coeff_dir: str = ""
    thresholds_file: str = ""


@dataclass
class AlgorithmConfig:
    """Algorithm toggle switches (1=enabled, 0=disabled)."""
    cloudmask: int = 1
    cloudamount: int = 0
    cloudphase: int = 0
    cloudtopz: int = 0
    cloudtau_day: int = 0
    cloudtau_night: int = 0
    cloudtype_ii: int = 0
    surface_sst: int = 0
    write_intermediate: int = 0


@dataclass
class NWPConfig:
    """NWP data source configuration."""
    source_id: int = 5
    rtm_option: int = 0
    grib_file_1: str = ""
    grib_file_2: str = ""


@dataclass
class ProcessingConfig:
    """Batch processing parameters."""
    year_start: int = 2020
    month_start: int = 1
    day_start: int = 1
    year_end: int = 2020
    month_end: int = 1
    day_end: int = 1
    hour_min: int = 0
    hour_max: int = 24
    n_threads: int = 1


@dataclass
class FY3Config:
    """Top-level configuration container."""
    sensor: SensorConfig = field(default_factory=SensorConfig)
    paths: PathConfig = field(default_factory=PathConfig)
    algorithms: AlgorithmConfig = field(default_factory=AlgorithmConfig)
    nwp: NWPConfig = field(default_factory=NWPConfig)
    processing: ProcessingConfig = field(default_factory=ProcessingConfig)

    # Input data files (for single-orbit processing)
    geo_file: str = ""
    l1b_file: str = ""
    oisst_file: str = ""

    # Output product files
    output_clm: str = ""
    output_cla: str = ""
    output_clp: str = ""
    output_ctp: str = ""
    output_cot: str = ""
    output_con: str = ""
    output_sst: str = ""
    output_intermediate: str = ""

    def to_namelist(self, filepath: str) -> None:
        """Write config as Fortran namelist format for backward compatibility."""
        lines = [
            "&config",
            f'  fylat_sensor_id     = {self.sensor.sensor_id},',
            f'  code_root_path      = "{self.paths.code_root}",',
            f'  L1b_data_path       = "{self.paths.l1b_data_dir}",',
            f'  nwp_data_path       = "{self.paths.nwp_data_dir}",',
            f'  oisst_data_path     = "{self.paths.oisst_data_dir}",',
            f'  fy3_mersi_GEO_data  = "{self.geo_file}",',
            f'  fy3_mersi_L1b_data  = "{self.l1b_file}",',
            f'  fy3_mersi_CLM_data  = "{self.output_clm}",',
            f'  fy3_mersi_CLA_data  = "{self.output_cla}",',
            f'  fy3_mersi_CLP_data  = "{self.output_clp}",',
            f'  fy3_mersi_CTP_data  = "{self.output_ctp}",',
            f'  fy3_mersi_COT_data  = "{self.output_cot}",',
            f'  fy3_mersi_CON_data  = "{self.output_con}",',
            f'  fy3_mersi_SST_data  = "{self.output_sst}",',
            f'  fy3_intermediate    = "{self.output_intermediate}",',
            f'  fylat_nwp_opt       = {self.nwp.source_id},',
            f'  fylat_rtm_opt       = {self.nwp.rtm_option},',
            f'  nwp_grib_data1      = "{self.nwp.grib_file_1}",',
            f'  nwp_grib_data2      = "{self.nwp.grib_file_2}",',
            f'  oisst_data          = "{self.oisst_file}",',
            f'  cloudmask_id        = {self.algorithms.cloudmask},',
            f'  cloudamount_id      = {self.algorithms.cloudamount},',
            f'  cloudphase_id       = {self.algorithms.cloudphase},',
            f'  cloudtopz_id        = {self.algorithms.cloudtopz},',
            f'  cloudtau_day_id     = {self.algorithms.cloudtau_day},',
            f'  cloudtau_night_id   = {self.algorithms.cloudtau_night},',
            f'  cloudtypeII_id      = {self.algorithms.cloudtype_ii},',
            f'  surface_sst_id      = {self.algorithms.surface_sst},',
            f'  write_inter_id      = {self.algorithms.write_intermediate}/',
        ]
        with open(filepath, 'w') as f:
            f.write('\n'.join(lines) + '\n')


def _merge_dict(base: dict, override: dict) -> dict:
    """Recursively merge override into base."""
    result = base.copy()
    for k, v in override.items():
        if k in result and isinstance(result[k], dict) and isinstance(v, dict):
            result[k] = _merge_dict(result[k], v)
        else:
            result[k] = v
    return result


def load_config(config_path: str, overrides: Optional[dict] = None) -> FY3Config:
    """Load configuration from YAML file with optional overrides.

    Args:
        config_path: Path to YAML configuration file.
        overrides: Optional dict of overrides (e.g., from CLI arguments).

    Returns:
        Populated FY3Config instance.
    """
    with open(config_path, 'r') as f:
        raw = yaml.safe_load(f) or {}

    if overrides:
        raw = _merge_dict(raw, overrides)

    cfg = FY3Config()

    # Sensor
    if 'sensor' in raw:
        s = raw['sensor']
        cfg.sensor.sensor_id = s.get('id', 21)
        cfg.sensor.n_elem = s.get('n_elem', 2048)
        cfg.sensor.n_line = s.get('n_line', 2000)

    # Paths
    if 'paths' in raw:
        p = raw['paths']
        cfg.paths.code_root = p.get('code_root', '')
        cfg.paths.l1b_data_dir = p.get('l1b_data_dir', '')
        cfg.paths.geo_data_dir = p.get('geo_data_dir', '')
        cfg.paths.nwp_data_dir = p.get('nwp_data_dir', '')
        cfg.paths.oisst_data_dir = p.get('oisst_data_dir', '')
        cfg.paths.output_dir = p.get('output_dir', '')
        cfg.paths.coeff_dir = p.get('coeff_dir', '')
        cfg.paths.thresholds_file = p.get('thresholds_file', '')

    # Algorithms
    if 'algorithms' in raw:
        a = raw['algorithms']
        cfg.algorithms.cloudmask = a.get('cloudmask', 1)
        cfg.algorithms.cloudamount = a.get('cloudamount', 0)
        cfg.algorithms.cloudphase = a.get('cloudphase', 0)
        cfg.algorithms.cloudtopz = a.get('cloudtopz', 0)
        cfg.algorithms.cloudtau_day = a.get('cloudtau_day', 0)
        cfg.algorithms.cloudtau_night = a.get('cloudtau_night', 0)
        cfg.algorithms.cloudtype_ii = a.get('cloudtype_ii', 0)
        cfg.algorithms.surface_sst = a.get('surface_sst', 0)
        cfg.algorithms.write_intermediate = a.get('write_intermediate', 0)

    # NWP
    if 'nwp' in raw:
        n = raw['nwp']
        cfg.nwp.source_id = n.get('source_id', 5)
        cfg.nwp.rtm_option = n.get('rtm_option', 0)
        cfg.nwp.grib_file_1 = n.get('grib_file_1', '')
        cfg.nwp.grib_file_2 = n.get('grib_file_2', '')

    # Processing
    if 'processing' in raw:
        pr = raw['processing']
        cfg.processing.year_start = pr.get('year_start', 2020)
        cfg.processing.month_start = pr.get('month_start', 1)
        cfg.processing.day_start = pr.get('day_start', 1)
        cfg.processing.year_end = pr.get('year_end', 2020)
        cfg.processing.month_end = pr.get('month_end', 1)
        cfg.processing.day_end = pr.get('day_end', 1)
        cfg.processing.hour_min = pr.get('hour_min', 0)
        cfg.processing.hour_max = pr.get('hour_max', 24)
        cfg.processing.n_threads = pr.get('n_threads', 1)

    # Input files
    cfg.geo_file = raw.get('geo_file', '')
    cfg.l1b_file = raw.get('l1b_file', '')
    cfg.oisst_file = raw.get('oisst_file', '')

    # Output files
    cfg.output_clm = raw.get('output_clm', '')
    cfg.output_cla = raw.get('output_cla', '')
    cfg.output_clp = raw.get('output_clp', '')
    cfg.output_ctp = raw.get('output_ctp', '')
    cfg.output_cot = raw.get('output_cot', '')
    cfg.output_con = raw.get('output_con', '')
    cfg.output_sst = raw.get('output_sst', '')
    cfg.output_intermediate = raw.get('output_intermediate', '')

    return cfg
