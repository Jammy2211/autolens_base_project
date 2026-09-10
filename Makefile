.PHONY: validate-literature-citations audit test

validate-literature-citations:
	python -m autoassistant.literature validate-citations

# Symbol/idiom audit of skills + wiki against the installed stack.
audit:
	python autoassistant/audit_skill_apis.py

# Assistant tooling test suite (slow: the gate tests import autolens per case).
test:
	python -m pytest autoassistant/tests -q
