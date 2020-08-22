def create_configuration(name: str = 'LOWBD2', **kwargs):
    from muser.data_models.parameters import muser_path
    import os
    conf_dir = muser_path('configurations')

    location = EarthLocation(lon=115.2505 * u.deg, lat=42.211833333 * u.deg, height=1365.0 * u.m)
    if name == 'MUSER2':
        antfile = os.path.join(conf_dir, 'muser-2.csv')
        lowcore = create_configuration_from_file(antfile=antfile,
                                                 mount='altaz', names='MUSER_%d',
                                                 diameter=2.0, name='MUSER', location=location, **kwargs)
    else:
        antfile = os.path.join(conf_dir, 'muser-1.csv')
        lowcore = create_configuration_from_file(antfile=antfile,
                                                 mount='altaz', names='MUSER_%d',
                                                 diameter=4.5, name='MUSER', location=location, **kwargs)
    return lowcore