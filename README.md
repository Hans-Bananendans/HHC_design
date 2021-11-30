# HHC_design
A codebase for the design and modelling of Helmholtz cages.

hhc_classes.py contains a Helmholtz cage class **HHCage** which can be filled up with three **HHCoil** objects. Each **HHCoil** object needs a **PowerSupply** object to properly evaluate the performance.

Quick instructions:
 - See *hhc_tests.py* for an overview of how the classes are constructed, and the methods are called. 
 - Each class implemented in *hhc_classes.py* contains a method called *listmethods()* which can be used to list all the methods that the class has (except those prefixed by "__").
