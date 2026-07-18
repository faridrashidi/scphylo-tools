"""Verify package logging behavior."""

import pytest

import scphylo as scp


class TestLogging:
    """Exercise all public logging levels."""

    def test_logging(self):
        """Verify normal messages and error escalation across log levels."""
        scp.logg.print_version()
        scp.logg.debug("DEBUG")
        scp.logg.hint("HINT")
        scp.logg.info("INFO")
        scp.logg.warn("WARN")
        scp.logg.info("TIME", time=True, color="red")
        with pytest.raises(RuntimeError):
            scp.logg.error("ERROR")
        assert True
