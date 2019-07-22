from pygslib.surpac import Surpac

def test_surpac():    
    SURPAC = Surpac('C:/OG_Python/pygslib/tests/test_surpac/surpac_strings.txt')
    assert(SURPAC.location == 'smp')
    assert(SURPAC.date == '8-Feb-94')
    assert(SURPAC.purpose == '')
    assert(SURPAC.memo == 'Sample string file')
    assert(SURPAC.axis.x_1 == 2440.000)
    assert(SURPAC.axis.z_2 == 152.541)
    assert(SURPAC.records[0].string_record[5].d[1] == 'second description field')
    assert(SURPAC.records[1].string_record[8].x == 2449.300)
    assert(len(SURPAC.records[2].string_record) == 15)
    assert(len(SURPAC.records) == 4)
    assert(SURPAC.records[3].number == 93 )
    assert(SURPAC.records[0].string_record[12].length_d == 1)
    # assert(SURPAC.records[0].string_record[13].length_d == 0)
    assert(SURPAC.records[0].string_record[5].length_d == 3)
test_surpac()