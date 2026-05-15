#!/usr/bin/env python3
"""Convert Fortran ASCII threshold files to YAML format.

Usage:
    python convert_thresholds.py <input_file> <output_file>

Example:
    python convert_thresholds.py ../coeff/fylat_thresholds.mersi.ii3d.v8 \
        ../config/thresholds/mersi_ii3d_v8.yaml
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import yaml


def parse_threshold_file(filepath: str) -> dict[str, list[float]]:
    """Parse a Fortran threshold file with 'NAME : VALUE' format.

    Lines starting with '!' are comments. Values can be comma-separated.
    Name matching is case-insensitive and whitespace-insensitive.
    """
    thresholds = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip empty lines and comments
            if not line or line.startswith('!'):
                continue
            # Skip non-data lines (like rcs_id, version strings)
            if ':' not in line:
                continue
            # Split on first colon only
            parts = line.split(':', 1)
            if len(parts) != 2:
                continue
            name = parts[0].strip().lower()
            value_str = parts[1].strip()
            # Strip inline comments (e.g., "1.0  ! comment")
            if '!' in value_str:
                value_str = value_str[:value_str.index('!')].strip()
            # Skip string values
            if name in ('rcs_id', 'thresholds_file_ver'):
                continue
            # Parse numeric values (comma-separated)
            values = []
            for v in value_str.split(','):
                v = v.strip()
                if v:
                    try:
                        values.append(float(v))
                    except ValueError:
                        # Skip non-numeric values
                        pass
            if values:
                thresholds[name] = values
    return thresholds


def organize_thresholds(raw: dict[str, list[float]]) -> dict:
    """Organize flat threshold dict into structured categories."""

    def get(key: str, default=None):
        return raw.get(key, default)

    def get4(key: str):
        """Get 4-element threshold array [locut, midpt, hicut, power]."""
        v = raw.get(key)
        if v and len(v) >= 4:
            return [v[0], v[1], v[2], v[3]]
        return v

    def get1(key: str):
        """Get single-element threshold."""
        v = raw.get(key)
        if v and len(v) >= 1:
            return v[0]
        return None

    def get2(key: str):
        """Get 2-element threshold array."""
        v = raw.get(key)
        if v and len(v) >= 2:
            return [v[0], v[1]]
        return v

    def get3(key: str):
        """Get 3-element threshold array."""
        v = raw.get(key)
        if v and len(v) >= 3:
            return [v[0], v[1], v[2]]
        return v

    organized = {
        "snow_mask": {
            "bt11_threshold": get1("sm_bt11"),
            "ndsi_threshold": get1("sm_ndsi"),
            "ref086_threshold": get1("sm_ref2"),
            "ref138_threshold": get1("sm_ref3"),
            "co2_threshold": get1("sm_co2"),
            "btd_85_11_threshold": get1("sm85_11"),
            "btd_37_11_threshold": get1("sm37_11"),
            "btd_37_11_hel_threshold": get1("sm37_11hel"),
            "nir_threshold": get1("sm_mnir"),
        },
        "land_day": {
            "vrat": get4("dlvrat"),
            "ref064": get4("dlref1"),
            "ref138": get4("dlref3"),
            "btd_11_4": get4("dl11_4lo"),
            "btd_11_12": get1("dl11_12hi"),
            "co2": get4("dlco2"),
            "h2o": get4("dlh20"),
            "tci": get2("dltci"),
        },
        "land_day_coast": {
            "ref064": get4("dlref1_t2"),
            "ref138": get4("dlref3_t2"),
            "btd_11_4": get4("dl11_4lo_t2"),
            "btd_11_12": get1("dl11_12hi_t2"),
            "co2": get4("dlco2_t2"),
            "h2o": get4("dlh20_t2"),
            "tci": get2("dltci_t2"),
        },
        "land_day_desert": {
            "btd_11_12": get4("lds11_12hi"),
            "ref086": get4("ldsref2"),
            "ref138": get4("ldsref3"),
            "btd_11_4_lo": get4("lds11_4lo"),
            "btd_11_4_hi": get4("lds11_4hi"),
            "co2": get4("ldsco2"),
            "h2o": get4("ldsh20"),
            "tci": get2("ldstci"),
        },
        "land_day_desert_coast": {
            "btd_11_12": get4("lds11_12hi_c"),
            "btd_11_4_hi": get4("lds11_4hi_c"),
            "btd_11_4_lo": get4("lds11_4lo_c"),
            "co2": get4("ldsco2_c"),
            "h2o": get4("ldsh20_c"),
            "ref086": get4("ldsref2_c"),
            "ref138": get4("ldsref3_c"),
            "tci": get2("ldstci_c"),
        },
        "land_day_polar": {
            "btd_11_12": get1("pdl11_12hi"),
            "btd_11_4": get4("pdl11_4lo"),
            "h2o": get4("pdlh20"),
            "ref064": get4("pdlref1"),
            "ref138": get4("pdlref3"),
            "vrat": get4("pdlvrat"),
            "tci": get2("pdltci"),
        },
        "land_day_polar_coast": {
            "btd_11_12": get1("pdl11_12hi_t2"),
            "btd_11_4": get4("pdl11_4lo_t2"),
            "h2o": get4("pdlh20_t2"),
            "ref064": get4("pdlref1_t2"),
            "ref138": get4("pdlref3_t2"),
            "tci": get2("pdltci_t2"),
        },
        "land_nite": {
            "btd_38_12_hi": get4("nl4_12hi"),
            "btd_38_12_lo": get4("nl4_12lo"),
            "co2": get4("nlco2"),
            "h2o": get4("nlh20"),
            "btd_73_11": get4("nl7_11s"),
            "btd_11_38_lo": get4("nl_11_4l"),
            "btd_11_38_hi": get4("nl_11_4h"),
            "btd_11_38_mid": get4("nl_11_4m"),
            "bt_diff_bounds": get2("bt_diff_bounds"),
            "btd_11_12": get1("nl11_12hi"),
        },
        "land_nite_polar": {
            "h2o": get4("pnlh20"),
            "btd_11_12": get1("pnl11_12hi"),
        },
        "ocean_day": {
            "btd_11_12": get1("do11_12hi"),
            "btd_11_4": get4("do11_4lo"),
            "bt11": get4("dobt11"),
            "co2": get4("doco2"),
            "h2o": get4("doh20"),
            "ref086": get4("doref2"),
            "ref138": get4("doref3"),
            "vrat_hi": get4("dovrathi"),
            "vrat_lo": get4("dovratlo"),
            "tci": get2("dotci"),
        },
        "ocean_day_polar": {
            "btd_11_12": get1("pdo11_12hi"),
            "btd_11_4": get4("pdo11_4lo"),
            "bt11": get4("pdobt11"),
            "h2o": get4("pdoh20"),
            "ref086": get4("pdoref2"),
            "ref138": get4("pdoref3"),
            "vrat_hi": get4("pdovrathi"),
            "vrat_lo": get4("pdovratlo"),
            "tci": get2("pdotci"),
        },
        "ocean_day_spatial_var": {
            "var_11um": get1("dovar11"),
        },
        "ocean_nite": {
            "btd_11_12": get1("no11_12hi"),
            "btd_11_4": get4("no11_4lo"),
            "bt11": get4("nobt11"),
            "co2": get4("noco2"),
            "h2o": get4("noh20"),
            "btd_86_73": get4("no86_73"),
            "var_11um": get4("no_11var"),
        },
        "ocean_nite_polar": {
            "btd_11_12": get1("pno11_12hi"),
            "btd_11_4": get4("pno11_4lo"),
            "bt11": get4("pnobt11"),
            "h2o": get4("pnoh20"),
            "btd_86_73": get4("pno86_73"),
            "var_11um": get4("pno_11var"),
        },
        "day_snow": {
            "btd_11_12": get1("ds11_12hi"),
            "btd_38_11": get4("ds4_11"),
            "btd_38_11_hel": get4("ds4_11hel"),
            "co2": get4("dsco2"),
            "h2o": get4("dsh20"),
            "ref138": get4("dsref3"),
            "tci": get2("dstci"),
            "btd_11_12_adj": get1("ds11_12adj"),
        },
        "day_snow_polar": {
            "btd_11_12": get1("dps11_12hi"),
            "h2o": get4("dpsh20"),
            "ref064": get4("dpsref1"),
            "ref138": get4("dpsref3"),
            "tci": get2("dpstci"),
            "btd_11_12_adj": get1("dps11_12adj"),
            "btd_38_11_lo": get4("dps4_11l"),
            "btd_38_11_hi": get4("dps4_11h"),
            "btd_38_11_mid1": get4("dps4_11m1"),
            "btd_38_11_mid2": get4("dps4_11m2"),
            "btd_38_11_mid3": get4("dps4_11m3"),
            "bt_11_bounds": raw.get("bt_11_bnds3"),
        },
        "nite_snow": {
            "btd_11_12": get1("ns11_12hi"),
            "btd_11_4": get4("ns11_4lo"),
            "btd_38_12": get4("ns4_12hi"),
            "co2": get4("nsco2"),
            "h2o": get4("nsh20"),
            "ref_065_11": get1("n65_11"),
            "btd_11_12_adj": get1("ns11_12adj"),
        },
        "nite_snow_polar": {
            "btd_11_12": get1("pns11_12hi"),
            "btd_38_12_lo": get4("pn_4_12l"),
            "btd_38_12_hi": get4("pn_4_12h"),
            "btd_38_12_mid1": get4("pn_4_12m1"),
            "btd_38_12_mid2": get4("pn_4_12m2"),
            "btd_38_12_mid3": get4("pn_4_12m3"),
            "btd_73_11_lo": get4("pn_7_11l"),
            "btd_73_11_hi": get4("pn_7_11h"),
            "btd_73_11_mid1": get4("pn_7_11m1"),
            "btd_73_11_mid2": get4("pn_7_11m2"),
            "btd_73_11_mid3": get4("pn_7_11m3"),
            "btd_73_11_lo_w": get4("pn_7_11lw"),
            "btd_73_11_hi_w": get4("pn_7_11hw"),
            "btd_73_11_mid1_w": get4("pn_7_11m1w"),
            "btd_73_11_mid2_w": get4("pn_7_11m2w"),
            "btd_73_11_mid3_w": get4("pn_7_11m3w"),
            "h2o": get4("pnsh20"),
            "btd_11_38_lo": get4("pn_11_4l"),
            "btd_11_38_hi": get4("pn_11_4h"),
            "btd_11_38_mid1": get4("pn_11_4m1"),
            "btd_11_38_mid2": get4("pn_11_4m2"),
            "btd_11_38_mid3": get4("pn_11_4m3"),
            "bt_11_bounds": raw.get("bt_11_bounds"),
            "bt_11_bounds2": raw.get("bt_11_bnds2"),
            "ref_065_11": get1("pn65_11"),
            "ref_138_11": get1("pn13_11"),
            "ref_073_11": get1("pn7_11"),
            "btd_11_12_adj": get1("pn11_12adj"),
        },
        "antarctic_day": {
            "h2o": get4("anth20"),
            "btd_38_11_lo": get4("ant4_11l"),
            "btd_38_11_hi": get4("ant4_11h"),
            "btd_38_11_mid1": get4("ant4_11m1"),
            "btd_38_11_mid2": get4("ant4_11m2"),
            "btd_38_11_mid3": get4("ant4_11m3"),
            "bt_11_bounds": raw.get("bt_11_bnds4"),
        },
        "pfmft": {
            "bt_11_max": get1("pfmft_11maxthre"),
            "btd_min": get1("pfmft_btd_min"),
            "land": get4("pfmft_land"),
            "ocean": get4("pfmft_ocean"),
            "snow": get4("pfmft_snow"),
            "cold": get4("pfmft_cold"),
        },
        "nfmft": {
            "max_threshold": get1("nfmft_maxthre"),
            "land": get4("nfmft_land"),
            "ocean": get4("nfmft_ocean"),
            "snow": get4("nfmft_snow"),
            "desert": get4("nfmft_desert"),
        },
        "land_restoral": {
            "r5_4_desert": get1("ldsr5_4_thr"),
            "r5_4_land": get1("ldr5_4_thr"),
            "btd_38_40": get1("ld20m22"),
            "btd_38_11": get1("ld22m31"),
            "bt11_hot_desert": get3("ldsbt11"),
            "bt11_hot_desert_bd": get3("ldsbt11bd"),
            "bt11_hot_nite": get3("lnbt11"),
        },
        "shadows": {
            "ref_nir": get2("shadnir"),
            "vrat": get1("shavrat"),
            "ref_124": get1("shad124"),
        },
        "sunglint": {
            "vrat": get2("snglntv"),
            "vrat_ch": get2("snglntvch"),
            "vrat_cl": get2("snglntvcl"),
            "btd_lo": get1("sg_tbdfl"),
            "btd_hi": get1("sg_tbdfh"),
            "ratio": get1("snglrat"),
            "ref_0deg": get4("snglnt0"),
            "ref_10deg": get4("snglnt10"),
            "ref_20deg": get4("snglnt20"),
            "bounds": get4("snglnt_bounds"),
        },
        "noncloud_obstruction": {
            "bt_38": get1("nc_bt37"),
            "btd_38_11": get1("nc37_11"),
            "ref_213": get1("nc21"),
            "btd_11_12": get1("nc11_12"),
            "ratio": get1("ncrat"),
            "vrat": get1("ncvrat"),
            "sigma": get1("ncsig"),
        },
        "ndvi_coast_swamp": {
            "bounds": get2("swc_ndvi"),
        },
    }

    return organized


def main():
    if len(sys.argv) < 3:
        print("Usage: python convert_thresholds.py <input_file> <output_file>")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    print(f"Reading thresholds from: {input_path}")
    raw = parse_threshold_file(input_path)
    print(f"  Found {len(raw)} threshold parameters")

    organized = organize_thresholds(raw)

    # Create output directory if needed
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        yaml.dump(organized, f, default_flow_style=False, sort_keys=False, allow_unicode=True)

    print(f"Written organized thresholds to: {output_path}")


if __name__ == '__main__':
    main()
