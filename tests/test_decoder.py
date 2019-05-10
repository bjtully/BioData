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


def test_command_line_default_viz(script_runner, tmp_path):
    p = tmp_path / "NORP_subset.txt"
    with open('tests/NORP_subset.txt', 'r') as f:
        p.write_text(f.read())

    ret = script_runner.run('KEGG-decoder',
                            '-i', str(p),
                            '-o', 'test.txt',
                            cwd=str(tmp_path))

    print(ret.stdout)
    print(ret.stderr)
    assert ret.success
    #assert ret.stdout == ''
    #assert ret.stderr == ''


def test_command_line_static(script_runner, tmp_path):
    p = tmp_path / "NORP_subset.txt"
    with open('tests/NORP_subset.txt', 'r') as f:
        p.write_text(f.read())

    ret = script_runner.run('KEGG-decoder',
                            '-i', str(p),
                            '-o', 'test.txt',
                            '-v', 'static',
                            cwd=str(tmp_path))

    print(ret.stdout)
    print(ret.stderr)
    assert ret.success
    #assert ret.stdout == ''
    #assert ret.stderr == ''


def test_command_line_interactive(script_runner, tmp_path):
    p = tmp_path / "NORP_subset.txt"
    with open('tests/NORP_subset.txt', 'r') as f:
        p.write_text(f.read())

    ret = script_runner.run('KEGG-decoder',
                            '-i', str(p),
                            '-o', 'test.txt',
                            '-v', 'interactive',
                            cwd=str(tmp_path))
    print(ret.stdout)
    print(ret.stderr)
    assert ret.success
    #assert ret.stdout == ''
    #assert ret.stderr == ''
