"""Fail closed if a committed SPOC source product changes."""

import hashlib
from pathlib import Path


EXPECTED = {
    "tess2019058134432-s0009-0000000035516889-0139-s_lc.fits": "70a40881930e7b374632a7afabb34ba2b18ff81841f6a622f03cf1d0268a06c0",
    "tess2021065132309-s0036-0000000035516889-0207-s_lc.fits": "6e465e081090b36cdeb0d467bc27270e99801ca3d3ee75b10946a741d9a5a1a1",
    "tess2023043185947-s0062-0000000035516889-0254-s_lc.fits": "c96bb25dadab40209390793713ad8ca6a579a0b224b06b9889f689cce555599a",
    "tess2023069172124-s0063-0000000035516889-0255-s_lc.fits": "8d61f6d1c3bc90f061ae4ce37ce9509ed52f4836022576ef6581510061a5601f",
    "tess2025042113628-s0089-0000000035516889-0286-s_lc.fits": "8f1e4763f361e85c0fe57036666264cc5a656cf2616e0fb1465497b3f1dd2e2f",
    "tess2025071122000-s0090-0000000035516889-0287-s_lc.fits": "c0e055be676307e1c19b9d9f51f9e879c739d210a421ec81c1e17ed63aed2428",
    "tess2026005125623-s0099-0000000035516889-0300-s_lc.fits": "64fd1ea3c5f66466b7b9c7400c40a96218ece4fd978f2627398617e5fec5cf90",
}


def test_complete_spoc_inventory_matches_frozen_checksums():
    data = Path(__file__).resolve().parents[1] / "data"
    observed = {path.name for path in data.glob("tess*_lc.fits")}
    assert observed == set(EXPECTED)
    for name, digest in EXPECTED.items():
        assert hashlib.sha256((data / name).read_bytes()).hexdigest() == digest
