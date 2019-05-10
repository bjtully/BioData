from KEGGDecoder import KEGG_decoder


def test_nitrogen():
    out = KEGG_decoder.nitrogen("K00368")
    assert 'nitrite reduction' in out
    assert out['nitrite reduction'] == 1

def test_arsenic():
    out = KEGG_decoder.arsenic(["K00537", "K03325"])
    assert 'Arsenic reduction' in out
    assert out['Arsenic reduction'] >= 0 
    assert out['Arsenic reduction'] <= 1 
    # there are four arsenic KOs
    # should be divisible by .25
    # this line relies on the implementation
    # which may not be ideal in the long run
    assert out['Arsenic reduction'] % .25 == 0   
