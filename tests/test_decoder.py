from KEGGDecoder import KEGG_decoder


def test_nitrogen():
    out = KEGG_decoder.nitrogen("K00368")
    assert 'nitrite reduction' in out
    assert out['nitrite reduction'] == 1
