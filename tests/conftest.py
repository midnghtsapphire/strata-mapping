"""Pytest configuration for strata-mapping."""


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "live_smoke: opt-in TNM and OpenTopography network tests; "
        "skipped unless STRATA_LIVE_SMOKE=1. OpenTopography also skips when "
        "OPENTOPO_API_KEY is unset. Never commit that key.",
    )
