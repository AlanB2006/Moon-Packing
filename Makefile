.PHONY: test figures verify clean

test:
	python -m pytest

figures:
	moon-packing-figures --figures 2 3 4 8 9 10

verify:
	sha256sum --check checksums.sha256

clean:
	rm -rf .pytest_cache build dist src/*.egg-info
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
