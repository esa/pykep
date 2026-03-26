.. _epoch:

Epoch class
################

The representation of an epoch, that is of a specific point in time, be it in the future or in the past, can be rather confusing. 
In `pykep` we opted to offer a dedicated class called :class:`~pykep.epoch` to offer a simple interface and, under the hoods, interfacing
seamlessly both to the c++ `std::chrono <https://en.cppreference.com/w/cpp/header/chrono>`_ 
library and to the python `datetime <https://docs.python.org/3/library/datetime.html>`_ module. 

.. note::
    In `pykep` the default Julian Date is the Modified Julian Date, defined as a `float` representing the number
    of days since the start of 2000-1-1. 

.. note::
    The date in `pykep` **does** account for leap seconds. If the user wishes to use the exact ISO 8601 representation of some epoch, 
    also including leap seconds, he will have to account for the offset himself. As of of 2023 this may account to maximum 28 seconds.
    `More info <https://en.wikipedia.org/wiki/Leap_second>`_ on leap seconds.

.. currentmodule:: pykep

-----------------------------------

.. autoclass:: pykep.epoch
   :members:
