
from pathlib import Path 
import pytest 

@pytest.fixture
def examples_dir(): 
    test_dir = Path(__file__).parent 
    return test_dir / 'example' 

@pytest.fixture 
def simple_fasta(examples_dir): 
    return examples_dir / 'simple.fasta'

@pytest.hookimpl()
def pytest_sessionfinish(session: pytest.Session, exitstatus):
    print(f"Collected {session.testscollected}, Failed {session.testsfailed}")
    for item in session.items: 
        print(item.reportinfo)
    print(f"Exit status: {exitstatus}")

    