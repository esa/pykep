from PyKEP import __extensions__

if (__extensions__['scipy']):
    from ._leg import leg
