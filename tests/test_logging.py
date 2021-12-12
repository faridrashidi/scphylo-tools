import pytest

import scphylo as scp


class TestLogging:
    def test_logging(self):
        scp.logg.print_version()
        scp.logg.debug("DEBUG")
        scp.logg.hint("HINT")
        scp.logg.info("INFO")
        scp.logg.warn("WARN")
        scp.logg.info("TIME", time=True, color="red")
        with pytest.raises(RuntimeError):
            scp.logg.error("ERROR")
        assert True
