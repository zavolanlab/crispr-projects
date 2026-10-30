from lincdisease import literature


def test_usable_alias():
    approved = {"TMEM99X"}
    assert not literature.usable_alias("ALIEN", approved)      # a word
    assert not literature.usable_alias("LEU2", approved)       # short, yeast gene
    assert not literature.usable_alias("TMEM99X", approved)    # another gene's symbol
    assert literature.usable_alias("LINC00657", approved)


def test_queries_are_exact_phrases():
    q = literature.build_queries("LINC01090")
    assert '(TITLE:"LINC01090" OR ABSTRACT:"LINC01090")' in q["n_liver"]
    assert "Retracted Publication" in q["n_retracted"]
