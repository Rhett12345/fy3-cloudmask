"""FY-3D MERSI-II 输入文件匹配工具。

给定日期或自动扫描，找到 L1B + GEO + NWP 三者都有的时次，
返回可直接用于 backend_compare_and_viz.py 的文件路径组合。

用法：
    # 作为模块导入
    from find_matched_files import find_matched_triplets, find_best_nwp

    # 命令行查询某天所有可用时次
    python find_matched_files.py --date 20230606

    # 命令行列出所有可用日期的第一个时次
    python find_matched_files.py --all
"""

from __future__ import annotations

import argparse
import os
import re
from datetime import datetime, timedelta
from pathlib import Path

# ============================================================
# 目录配置 — 按实际情况修改
# ============================================================
MERSI_ROOT = Path('/data/Data_yuq/mersi_test')   # 子目录结构: YYYYMMDD/FY3D_*.HDF
NWP_ROOT   = Path('/data/nwp')              # 子目录结构: YYYYMMDD/ORG/gfs0p25_41L_*

# NWP 文件名模板，时次为 00/03/06/09/12/15/18/21
NWP_PATTERN = 'fnl_{date}_{hh}_00'
NWP_HOURS   = [0, 3, 6, 9, 12, 15, 18, 21]   # UTC, 整点时次


# ============================================================
# 内部工具
# ============================================================

def _parse_mersi_time(filename: str) -> datetime | None:
    """从 L1B/GEO 文件名解析 UTC 时间。

    文件名格式：FY3D_MERSI_GBAL_L1_YYYYMMDD_HHMM_1000M_MS.HDF
    """
    m = re.search(r'_(\d{8})_(\d{4})_', filename)
    if not m:
        return None
    try:
        return datetime.strptime(m.group(1) + m.group(2), '%Y%m%d%H%M')
    except ValueError:
        return None


def _find_nwp(date_str: str, obs_hour: int) -> Path | None:
    """找与观测时间最近的 NWP 文件（向前取整到最近时次）。

    Args:
        date_str: 'YYYYMMDD'
        obs_hour: 观测 UTC 小时数（0-23）

    Returns:
        NWP 文件路径，找不到返回 None。
    """
    # 找不超过 obs_hour 的最近 NWP 时次
    best_hh = None
    for hh in sorted(NWP_HOURS, reverse=True):
        if hh <= obs_hour:
            best_hh = hh
            break
    if best_hh is None:
        # obs_hour 比最小时次还小（例如凌晨0点前），取前一天最后时次
        prev_date = (datetime.strptime(date_str, '%Y%m%d') - timedelta(days=1)).strftime('%Y%m%d')
        best_hh = max(NWP_HOURS)
        date_str = prev_date

    fname = NWP_PATTERN.format(date=date_str, hh=f'{best_hh:02d}')
    nwp_path = NWP_ROOT / date_str / fname
    return nwp_path if nwp_path.exists() else None


def find_matched_triplets(date_str: str) -> list[dict]:
    """查找某天所有 L1B + GEO + NWP 齐全的时次。

    Args:
        date_str: 'YYYYMMDD'

    Returns:
        列表，每个元素为字典：
        {
            'datetime': datetime,
            'l1b': Path,
            'geo': Path,
            'nwp': Path,
        }
        按时间排序。
    """
    mersi_dir = MERSI_ROOT / date_str
    if not mersi_dir.exists():
        return []

    # 收集所有 1000M L1B 文件
    l1b_files: dict[datetime, Path] = {}
    geo_files: dict[datetime, Path] = {}

    for f in mersi_dir.iterdir():
        if not f.suffix.upper() == '.HDF':
            continue
        t = _parse_mersi_time(f.name)
        if t is None:
            continue
        if '1000M_MS' in f.name:
            l1b_files[t] = f
        elif 'GEO1K_MS' in f.name:
            geo_files[t] = f

    results = []
    for t, l1b in sorted(l1b_files.items()):
        if t not in geo_files:
            continue  # 没有对应 GEO 文件
        nwp = _find_nwp(date_str, t.hour)
        if nwp is None:
            continue  # 没有可用 NWP
        results.append({
            'datetime': t,
            'l1b': l1b,
            'geo': geo_files[t],
            'nwp': nwp,
        })

    return results


def find_all_available_dates() -> list[str]:
    """列出 MERSI_ROOT 下所有有数据的日期（YYYYMMDD）。"""
    if not MERSI_ROOT.exists():
        return []
    dates = sorted(
        d.name for d in MERSI_ROOT.iterdir()
        if d.is_dir() and re.fullmatch(r'\d{8}', d.name)
    )
    return dates


def find_best_triplet(date_str: str) -> dict | None:
    """返回某天第一个可用的文件三元组，找不到返回 None。"""
    triplets = find_matched_triplets(date_str)
    return triplets[0] if triplets else None


# ============================================================
# 命令行入口
# ============================================================

def _print_triplet(t: dict) -> None:
    print(f"  时间: {t['datetime'].strftime('%Y-%m-%d %H:%M UTC')}")
    print(f"  L1B: {t['l1b']}")
    print(f"  GEO: {t['geo']}")
    print(f"  NWP: {t['nwp']}")
    print()


def main():
    parser = argparse.ArgumentParser(description='FY-3D MERSI-II 文件匹配工具')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--date', metavar='YYYYMMDD',
                       help='查询某天所有可用时次')
    group.add_argument('--all', action='store_true',
                       help='列出所有可用日期的第一个时次')
    args = parser.parse_args()

    if args.date:
        triplets = find_matched_triplets(args.date)
        if not triplets:
            print(f'[{args.date}] 没有找到完整的 L1B+GEO+NWP 组合')
            return
        print(f'[{args.date}] 找到 {len(triplets)} 个可用时次:\n')
        for t in triplets:
            _print_triplet(t)

        # 输出可直接粘贴到脚本的变量
        first = triplets[0]
        print('# 粘贴到 backend_compare_and_viz.py 的路径变量（第一个时次）:')
        print(f"L1B_FILE = '{first['l1b']}'")
        print(f"GEO_FILE = '{first['geo']}'")
        print(f"NWP_FILE = '{first['nwp']}'")

    elif args.all:
        dates = find_all_available_dates()
        print(f'共找到 {len(dates)} 个日期，以下只显示每天第一个可用时次:\n')
        found = 0
        for date_str in dates:
            t = find_best_triplet(date_str)
            if t:
                found += 1
                print(f'[{date_str}]')
                _print_triplet(t)
        print(f'共 {found} 天有完整数据。')


if __name__ == '__main__':
    main()
