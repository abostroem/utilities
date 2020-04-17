from utilities_az import connect_to_supernova

db, cursor = connect_to_supernova.get_cursor(config_dir='/Volumes/quince_external/lcophotometry/')  

def test_one_results():
    x = cursor.execute('SELECT * from targets where id=2') 
    assert x == 1

def test_multiple_results():
    x = cursor.execute('SELECT * from targets')
    assert x>1

def test_no_results():
    x = cursor.execute('SELECT * from targets where id=-1') 
    assert x == 0
