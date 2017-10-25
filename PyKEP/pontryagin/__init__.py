from PyKEP import __extensions__

if (__extensions__['scipy'] and __extensions__['mplot3d']):
    from ._leg import leg
else:
    pass
