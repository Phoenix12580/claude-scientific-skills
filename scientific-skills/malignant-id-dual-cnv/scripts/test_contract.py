from malignant_id.contract import suggest_columns, validate_contract


def test_validate_requires_fields():
    errors = validate_contract({"schema_version": "1.0.0"})
    assert any("missing required field" in e for e in errors)
