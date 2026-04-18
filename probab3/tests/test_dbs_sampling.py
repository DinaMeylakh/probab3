"""Assert basic usage of the CLI. Sample an DBS triple, then show it's statistics and plot it. Assert no errors."""

import pytest
import tempfile
import os

from click.testing import CliRunner
from probab3.__main__ import cli


@pytest.fixture
def dbs_sampler_out_file(request):
    """Return a temporary file name for the output of the mt sampler."""

    out_file = tempfile.NamedTemporaryFile(delete=False)
    out_file.close()

    def teardown():
        os.unlink(out_file.name)

    request.addfinalizer(teardown)

    return out_file.name
    

def assert_data_usage(mocker, runner, out_file_name):
    """Givn sampled mt data, assert showing stats and plotting works."""

    show_stats_result = runner.invoke(cli, ['show-stats' , '--input-path', out_file_name, '--verbose'])

    assert show_stats_result.exit_code == 0

    show_stats_result = runner.invoke(cli, ['show-stats' , '--input-path', out_file_name, '--verbose', '--merger-only'])

    assert show_stats_result.exit_code == 0

    mocker.patch("probab3.commands.plot.plotter.plt.show")
    plot_result = runner.invoke(cli, ['plot' , '--input-path', out_file_name, '--N', '0'])

    assert plot_result.exit_code == 0


def test_equal_mass_basic(mocker, dbs_sampler_out_file):
    """Sample and test the basic equal mass dynamical binary sequence."""

    runner = CliRunner()
    run_dbs_result = runner.invoke(cli, ['sample-dbs' , '--output_path', dbs_sampler_out_file, '--Mc', '1e6', '--rh', '1.0', 
                                '--rc', '0.1', '--N', '1', 'basic', '--m1', '20', '--m2', '20', 
                                '--m3', '20'])

    assert run_dbs_result.exit_code == 0

    assert_data_usage(mocker, runner, dbs_sampler_out_file)


def test_equal_mass(mocker, dbs_sampler_out_file):
    """Sample and test the equal mass dynamical binary sequence with stars to BHs pre evolution."""

    runner = CliRunner()
    run_dbs_result = runner.invoke(cli, ['sample-dbs' , '--output_path', dbs_sampler_out_file, '--Mc', '1e6', '--rh', '1.0', 
                                        '--rc', '0.1', '--N', '1', 'equal-mass'])

    assert run_dbs_result.exit_code == 0

    assert_data_usage(mocker, runner, dbs_sampler_out_file)



def test_unqual_mass_sample(mocker, dbs_sampler_out_file):
    """Sample and test the unequal mass dynamical binary sequence with stars to BHs pre evolution."""
    runner = CliRunner()
    run_dbs_result = runner.invoke(cli, ['sample-dbs' , '--output_path', dbs_sampler_out_file, '--Mc', '1e6', '--rh', '1.0', 
                                        '--rc', '0.1', '--N', '1', 'unequal-mass'])

    assert run_dbs_result.exit_code == 0

    assert_data_usage(mocker, runner, dbs_sampler_out_file)

