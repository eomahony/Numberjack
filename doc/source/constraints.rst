Constraints Reference
=====================


Single operand predicates
-------------------------

.. autoclass:: Numberjack.Predicate
   :members:
   :inherited-members:

.. autoclass:: Numberjack.Neg
.. autoclass:: Numberjack.Abs


Binary constraints
------------------

Each of the following classes take two sub-expressions or variables as
parameters and can usually be created by operator overloading. For example `x *
y` will create a :class:`Numberjack.Mul` instance if either `x` or `y` are
Numberjack variables.

.. autoclass:: Numberjack.BinPredicate
.. autoclass:: Numberjack.And
.. autoclass:: Numberjack.Or
.. autoclass:: Numberjack.Eq
.. autoclass:: Numberjack.Ne
.. autoclass:: Numberjack.Lt
.. autoclass:: Numberjack.Le
.. autoclass:: Numberjack.Gt
.. autoclass:: Numberjack.Ge
.. autoclass:: Numberjack.Mul
.. autoclass:: Numberjack.Div
.. autoclass:: Numberjack.Mod



Other constaints
----------------

.. autoclass:: Numberjack.Table



Optimisation
------------

.. autoclass:: Numberjack.Maximise
.. autofunction:: Numberjack.Maximize
.. autoclass:: Numberjack.Minimise
.. autofunction:: Numberjack.Minimize

