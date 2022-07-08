install: ## Install and check dependencies
	poetry build && pip3 install --force-reinstall dist/pgscatalog_utils-0.1.0-py3-none-any.whl