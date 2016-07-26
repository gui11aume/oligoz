#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Name: oligo.
Author: Guillaume Filion.
Date: October 10, 2012.
Version: 1.0.

MODULE SUMMARY:
The module contains 3 major sections.

The first section implements objects representing nucleic
acids and solutions of nucleic acids. This part models nucleic
acid annealing properties. The major class 'OligoSol' captures
the concepts of this section. An instance of 'OligoSol' has an
instance of 'Oligo' and hybridization conditions, which allows
to compute the Tm. Base-pairing logic is deferred to the oligo
which itself contains an instance of 'DNAseq'. The class
'DNAseq' (and 'ntSeq' in general) implement reverse comple-
mentation, whereas the thermodynamics property are implemented
in the class 'Oligo'.

The second section implements oligo QC checks. Conceptually,
oligo solutions need to be consistent and compatible. The first
asserts that oligos do not self-anneal and that they meet the
standards specified by the user (such as a minimum Tm). The
second asserts that pairs or sets of oligos must be compatible
in the same reaction (they must not anneal to each other for
instance).

The third section implements search engines. Those engines are
product-centric in the sense that they search sets of oligo
solutions that will amplify the given product in a PCR reaction.
The search itself is decomposed in two steps: find consistent
oligo solutions, from which sets of compatible oligo solutions
are then selected.

REFERENCES:
(1) Santalucia J. 1998. A unified view of polymer, dumbbell,
   and oligonucleotide DNA nearest-neighbor thermodynamics.
(2) Owczarzy et al. 2008. Predicting stability of DNA duplexes
   in solutions containing magnesium an monovalent cations.
   Biochemistry, 47:5336-53.

APPENDIX:
IUPAC degenerate DNA letters.
R = A+G
Y = C+T
M = C+A
K = G+T
W = A+T
S = C+G
B = C+G+T
V = C+G+A
D = G+A+T
H = C+A+T
N = A+G+C+T
"""

##############################################
##                                          ##
##                 IMPORT                   ##
##                                          ##
##############################################

import re
import errno
import sys
import warnings

from math import log, exp, sqrt
from string import maketrans
from collections import defaultdict


##############################################
##                                          ##
##                 GLOBALS                  ##
##                                          ##
##############################################

ABS0 = -273.15

STD_HYBCOND = {
      'conc': 1e-6,
      'mono': 50e-3,
      'di'  : 1.5e-3,
}

STD_RULES = {
      'min_length':        17,
      'max_length':        30,
      'min_Tm':   273.15 + 58,
      'max_Tm':   273.15 + 65,
      'max_self_anneal':    3,
      'max_delta_Tm':       3,
      'max_hetero_anneal':  3,
      'max_lanneal':        4,
}

qPCR_RULES = {
      'min_length':        20,
      'max_length':        24,
      'min_Tm':   273.15 + 58,
      'max_Tm':   273.15 + 62,
      'max_self_anneal':    2,
      'max_delta_Tm':       4,
      'max_hetero_anneal':  3,
      'max_lanneal':        4,
}

RELAXED_RULES = {
      'min_length':        15,
      'max_length':        35,
      'min_Tm':   273.15 + 48,
      'max_Tm':   273.15 + 75,
      'max_self_anneal':    6,
      'max_delta_Tm':      15,
      'max_hetero_anneal':  8,
      'max_lanneal':       20,
}


##############################################
##                                          ##
##        EXCEPTIONS AND WARNINGS           ##
##                                          ##
##############################################


class ForbiddenCharacter(Exception):
   pass

def _warn(msg, category=UserWarning, filename='', lineno=-1):
   sys.stderr.write(str(msg))

# Monkey patch 'warnings.showwarning()'.
warnings.showwarning = _warn


##############################################
##                                          ##
##    CLASSES FOR DNA SEQUENCE OBJECTS      ##
##                                          ##
##############################################


##############################################
##########           ntSeq          ##########
##############################################

class ntSeq(str):
   """
   Abstract nucleic acid sequence.

   Abstract super class that provides the 'revcomp' property for
   reverse complementation. Extensions need to specify a 'bp'
   attribute and a 'validateSeq()' method.
   """

   ## DELEGATE ##
   def validateSeq(self, seq):
      raise NotImplementedError()

   # Overwrite __new__ because 'str' is immutable.
   def __new__(self, seq):
      seq = seq.upper()
      return str.__new__(self, seq)

   def __init__(self, seq):
      self.validateSeq(seq)


   # Avoid the confusion with instances of 'str'.
   # Use 'str(...)' to display the sequence.
   def __repr__(self):
      return '<%s.%s instance at %s>' % \
         (self.__module__, self.__class__.__name__, hex(id(self)))

   # TODO: The 'RevComp()' method is obsolete. Remove it after
   # checking that the module is stable.
   #def RevComp(self):
   #   """Return the Watson-Crick reverse complement of the
   #   sequence as a string.
   #   Requires a 'bp' attribute."""
   #
   #   return self.translate(self.bp)[::-1]

   @property
   def revcomp(self):
      """Return the Watson-Crick reverse complement of the
      sequence as a string.
      Requires a 'bp' attribute."""

      return self.translate(self.bp)[::-1]


##############################################
##########          DNAseq          ##########
##############################################

class DNAseq(ntSeq):
   """
   Non degenerate DNA sequence.
   Has a 'bp' attribute and a 'toRNA()' method.
   Inherits the 'revcomp' property from ntSeq.
   """

   # Complementary base-pairing.
   bp = maketrans('gatcGATC', 'ctagCTAG')

   def validateSeq(self, seq):
      # Check for non DNA letters.
      if self.__class__.isDNA(seq) is False:
         raise ForbiddenCharacter('non DNA letter(s) in ' + seq)

   def toRNA(self):
      return RNAseq(self.replace('T', 'U'))

   @classmethod
   def isDNA(self, seq):
      """Class method to check whether an object with a 'str' method
      consists only of DNA letters."""
      return re.search('[^GATC]', str(seq).upper()) is None


##############################################
##########          RNAseq          ##########
##############################################

class RNAseq(ntSeq):
   '''
   Non degenerate RNA sequence.
   Has a 'bp' attribute and a 'toDNA()' method.
   Inherits the 'revcomp' property from ntSeq.
   '''

   # Complementary base-pairing.
   bp = maketrans('gaucGAUC', 'cuagCUAG')

   def validateSeq(self, seq):
      # Check for non RNA letters.
      if re.search('[^GAUC]', seq):
         raise ForbiddenCharacter('non RNA letter(s) in ' + seq)


##############################################
##########         DNAdseq          ##########
##############################################

class DNAdseq(DNAseq):
   '''
   Degenerate DNA sequence.
   Inherits the 'toRNA()' method from DNAseq.
   '''
   # Complementary base-pairing.
   bp = maketrans('gatcrymkbvdhGATCRYMKBVDH', 'ctagyrkmvbhdCTAGYRKMVBHD')

   def validateSeq(self, seq):
      # Check for non degenerate DNA letters.
      if re.search('[^ATGCRYMKWSBVDHN]', seq):
         raise ForbiddenCharacter('non degenerate DNA letter(s) in ' + seq)

##############################################
##########         RNAdseq          ##########
##############################################

class RNAdseq(RNAseq):
   '''
   Degenerate RNA sequence.
   Inherits the 'toDNA()' method from RNAseq.
   '''

   # Complementary base-pairing.
   bp = maketrans('gaucrymkbvdhGAUCRYMKBVDH', 'cuagyrkmvbhdCUAGYRKMVBHD')

   def validateSeq(self, seq):
      # Check for non degenerate RNA letters.
      if re.search('[^AUGCRYMKWSBVDHN]', seq):
         raise ForbiddenCharacter('non degenerate RNA letter(s) in ' + seq)



##############################################
###########         Oligo          ###########
##############################################

class Oligo:
   """
   DNA oligonucleotide.

   The class attributes of the Oligo are physical constants for
   computing annealing properties of DNA strands. The 'seq'
   attribute is the part of the sequence that contributes to
   specific annealing, the 'extraNt' attribute is that part in 3'
   that does not contribute.

   The Tm is not an intrinsic property of an oligo because it
   depends on the chemical composition of the oligo solution.
   For this reason, the 'Tm' attribute is defined only for objects
   of class 'OligoSol'.

   Instance attributes:
     'seq': an instance of DNAseq.
     'extraNt': a 'str'.

   Instance methods:
     'misprime()': returns extensible 3' annealing on a template.

   Has-a interface requirements:
     'seq': an instance of 'DNAseq' (provided as a 'str' upon
       initialization) with a 'str()', '__getitem__()' and
       a 'revcomp' property.
     'extraNt': a 'str' containing DNA letters or N.
   """

   # Nearest-neighbor constants from ref (1).
   nnConst = {
   # H is in kcal/mol, S in cal/K mol.
   # Initiation terms.
    'A' : { 'H' : 2.3, 'S' : 4.1, },
    'C' : { 'H' : 0.1, 'S' : -2.8, },
    'G' : { 'H' : 0.1, 'S' : -2.8, },
    'T' : { 'H' : 2.3, 'S' : 4.1, },
   # Neares neighbor terms.
    'AA' : { 'H' : -7.9, 'S' : -22.2, },
    'AC' : { 'H' : -8.4, 'S' : -22.4, },
    'AG' : { 'H' : -7.8, 'S' : -21.0, },
    'AT' : { 'H' : -7.2, 'S' : -20.4, },
    'CA' : { 'H' : -8.5, 'S' : -22.7, },
    'CC' : { 'H' : -8.0, 'S' : -19.9, },
    'CG' : { 'H' : -10.6, 'S' : -27.2, },
    'CT' : { 'H' : -7.8, 'S' : -21.0, },
    'GA' : { 'H' : -8.2, 'S' : -22.2, },
    'GC' : { 'H' : -9.8, 'S' : -24.4, },
    'GG' : { 'H' : -8.0, 'S' : -19.9, },
    'GT' : { 'H' : -8.4, 'S' : -22.4, },
    'TA' : { 'H' : -7.2, 'S' : -21.3, },
    'TC' : { 'H' : -8.2, 'S' : -22.2, },
    'TG' : { 'H' : -8.5, 'S' : -22.7, },
    'TT' : { 'H' : -7.9, 'S' : -22.2, }
   }

   def __init__(self, seq, extraNt=''):
      """
      Initialize seq, and extraNt.

      Attribute 'thermoParams' is defined JIT upon first
      request.

      seq: an instance of 'DNAseq'.
      extraNt: a 'str' with DNA letters and free nucleotides
        indicated as a '?'.
      thermoParams: a dictionary containing the standard
        (i.e. in 1M NaCl) thermodynamic properties of the
        oligo:
        H: enthalpy of annealing in cal/mol.
        S: entropy of annealing in cal/K/mol.
      """

      self.seq = DNAseq(str(seq))
      self.extraNt = extraNt

   def __str__(self):
      return str(self.extraNt.lower() + self.seq)

   def __len__(self):
      return len(str(self))

   def __getattr__(self, name):
      """
      The thermodynamic attribute 'thermoParams' is defined
      just in time upon first request (seemless for the user).

      Compute enthalpy (H) and entropy (S) at a concentration
      of 1 M NaCl, according to the method presented in ref (1).
      H is in cal/mol and S in cal/K mol.
      """

      if name == 'thermoParams':
         H = Oligo.nnConst[self.seq[0]]['H'] + \
               Oligo.nnConst[self.seq[-1:]]['H']
         S = Oligo.nnConst[self.seq[0]]['S'] + \
               Oligo.nnConst[self.seq[-1:]]['S']
         # Iterate over dinucleotides.
         for i in range(len(self.seq)-1):
            H += Oligo.nnConst[self.seq[i:i+2]]['H']
            S += Oligo.nnConst[self.seq[i:i+2]]['S']
         self.thermoParams = {'H': 1000 * H, 'S': S}
         return self.thermoParams
      else:
         raise AttributeError()


   def misprime(self, seq=None, strict=True):
      """
      Check for mispriming between 'Oligo' instance and provided
      sequence (or self).

      RETURN:
         A 'str' with the longest 3' annealing of the 'Oligo'
         instance on the sequence defined by seq.

      ARGUMENTS:
         'seq': an object with a str() method, specifying a
           template for mispriming. Defaults to None, in which
           case the oligo itself is used for template.
         'strict': True or False, specifies whether only extendable
           mispriming are returned. Defaults to True.
      """

      if seq is None: seq = self.seq

      # Make sure that 'seq' is DNA.
      if DNAseq.isDNA(seq) is False:
         raise Exception("non DNA letters in '%s'" % seq)

      # Remove first nucleotide of sequence if 'strict' is set
      # to 'True'. In that case, mispriming is sure to be
      # extendable by at least one nucleotide. Otherwise, the
      # annealing can happen as shown in the example below,
      # whereby no extension is possible.
      #
      #               'Oligo'      5' TGCTGATGCGA 3'
      #               'seq'              5' ACGCT 5'
      #
      if strict is False:
         seq = str(seq).upper()
      else:
         seq = str(seq).upper()[1:]

      # This longest substring search is inefficient. Do not use
      # for long DNA sequences.
      rc = DNAseq(self.extraNt + self.seq).revcomp
      for i in range(1, len(rc)+1):
         if rc[:i] not in seq: break
      return rc[:i-1]


   def lanneal(self, seq=None, strict=True):
      """
      Compute the length of the maximum annealing between
      'Oligo' instance and the provided sequence (or self).

      RETURN:
         A 'str' with the longest annealing of the 'Oligo'
         instance on the sequence defined by seq.

      ARGUMENTS:
         'seq': an object with a str() method, specifying a
           template for mispriming. Defaults to None, in which
           case the oligo itself is used for template.
      """

      if seq is None: seq = self

      # Make sure that the sequence 'seq' is DNA.
      if DNAseq.isDNA(seq) is False:
         raise Exception("non DNA letters in '%s'" % seq)

      rc = DNAseq(self.extraNt + self.seq).revcomp

      # Find the longest common substring
      # by dynamic programming.
      lann = 0
      old_array = [0] * (len(rc)+1)
      for i, x in enumerate(str(seq).upper()):
         new_array = [0] * (len(rc)+1)
         for j, y in enumerate(rc):
            new_array[j+1] = old_array[j] + 1 if x == y else 0
         lann = max(lann, max(new_array))
         old_array = new_array

      return lann


##############################################
###########        OligoSol        ###########
##############################################

class OligoSol:
   """
   Oligo solution.

   Instance attributes:
   hybCond: a dictionary with hybridization conditions.
   oligo: an instance of Oligo.
   Tm0: the standard Nearest Neighbor Tm (in K) in 1 M NaCl.
   Tm: the predicted Tm (in K) in the given hybridization
     conditions.

   Instance methods:
   inverseTm(): computes the fraction of bound target at the
     given temperature.

   Has-a interface requirements:
   oligo: an instance of Oligo, with attributes 'seq' and
   hybCond: a dictionary containing:
     conc: oligo concentration in mol/L.
     mono: concentration of free monovalent cations in mol/L.
     di: concentration of free divalent (Mg2+) ions in mol/L.
       Note that dNTPs have a chelating effect on Mg2+ ions,
       so that [free Mg2+] = [Mg2+] - [dNTPs].
     target: (optional) target DNA concentration in mol/L
       (default: negligible)
    'thermoParams' .
   """

   # Owczarzy salt correction parameters from ref (2).
   saltP = {
      'a' : 3.92e-5,
      'b' : -9.11e-6,
      'c' : 6.26e-5,
      'd' : 1.42e-5,
      'e' : -4.82e-4,
      'f' : 5.25e-4,
      'g' : 8.31e-5
   }

   def __init__(self, oligo, hybCond=STD_HYBCOND):
      """Initialize the 'oligo' and 'hybCond' attributes.
      Do not initialize the 'Tm' attribute, which is defined
      just in time upon first request.

      oligo: instance of Oligo.
      hybCond: dictionary with keys:
         'conc': oligo concentration (mol/L).
         'mono': concentration of monovalent ions (mol/L).
         'di':   concentration of divalent cations (mol/L).
         'target': target concentration (mol/L), defaults to 0.

      Computing the Tm is somewhat expensive, so this allows to
      save time if the user does not actually need it."""

      # Keys 'conc', 'mono' and 'di' are required.
      # The key 'target' defaults to 0.
      missing = set(['conc', 'mono', 'di']).difference(hybCond)
      if missing:
         raise Exception('miss key(s): ' + ', '.join(missing))

      self.hybCond = hybCond
      self.hybCond.setdefault('target', 0)
      self.oligo = oligo if isinstance(oligo, Oligo) else Oligo(oligo)

   def __getattr__(self, name):
      """The attributes 'Tm' and 'Tm0' are defined just in time
      upon first request (seemless for the user)."""

      if name == 'Tm0':
         # Compute the nearest-neighbor Tm in K for [Na] = 1 mol/L.
         # By definition, Tm is the temperature for which half of
         # the target is annealed to the oligo at equilibrium. This
         # assumes [oligo] > [target].
         H = self.oligo.thermoParams['H']
         S = self.oligo.thermoParams['S']
         conc = self.hybCond['conc']
         target = self.hybCond['target']

         if target > conc:
            raise Exception ('Tm not defined given hybCond')

         # Tm = delta_H / (delta_S + R * log([A] - [B]/2))
         # For a palindrome, Tm = delta_H / (delta_S + R * log([A]))
         self.Tm0 = H / (S + 1.987 * log (conc - target/2))

         return self.Tm0

      elif name == 'Tm':
         # Compute the predicted Tm in K given hybridization
         # conditions. Use an additive correction term for 1/Tm.
         correction = self.OwczarzyCorrection()
         self.Tm = 1 / (1/self.Tm0 + correction)
         return self.Tm

      else:
         raise AttributeError(name)

   def __str__(self):
      return str(self.oligo)

   def __len__(self):
      return len(self.oligo)

   def OwczarzyCorrection (self):
      """Return an additive correction factor for 1/Tm.
      Use methods and notations described in ref (2)."""

      size = len(self.oligo.seq)
      f_GC = len(re.findall('G|C', self.oligo.seq)) / size
      mono = self.hybCond['mono']
      di = self.hybCond['di']

      if mono == 0 or sqrt(di) > 6 * mono:
         # R > 6: use formula (16) from ref (2).
         return OligoSol.saltP['a'] + OligoSol.saltP['b'] * log(di) +    \
            f_GC * (OligoSol.saltP['c'] + OligoSol.saltP['d']*log(di)) + \
            (OligoSol.saltP['e'] + OligoSol.saltP['f']*log(di) +         \
               OligoSol.saltP['g']*log(di))**2 / (2*(size- 1))

      if sqrt(di) > 0.22 * mono:
         # 0.22 < R < 6.0: use formula (16) and (18-20) from ref (2).
         updated_a = 3.92e-5 * (0.843 - 0.352 * sqrt(mono) * log(mono))
         updated_d = 1.42e-5 * (1.279 - 4.03e-3 * log(mono) -
            8.03e-3 * log(mono)**2);
         updated_g = 8.31e-5 * (0.486 - 0.258 * log(mono) +
            5.25e-3 * log(mono)**3);

         return updated_a + OligoSol.saltP['b'] * log(di) +       \
            f_GC * (OligoSol.saltP['c'] + updated_d*log(di)) +    \
            (OligoSol.saltP['e'] + OligoSol.saltP['f']*log(di) +  \
               updated_g*log(di))**2 / (2*(size- 1))

      # R < 0.22: use formula (4) from ref(2).
      return (4.29*f_GC-3.95) * 1e-5 * log(mono) + \
            9.40e-6 * (log(mono))**2

   def inverseTm(self, temp, Kelvin=False):
      """Return the fraction of bound target in the specified
      hybridization conditions."""

      if Kelvin is False:
         temp += 273.15

      # Recompute the entropy in the hybridization conditions.
      # Assume that salts do not affect enthalpy
      H0 = self.oligo.thermoParams['H']
      S0 = self.oligo.thermoParams['S']
      S = (1/self.Tm - 1/self.Tm0) * H0 + S0

      # Use formula k = exp(-dG/RT).
      k = exp((H0-temp*S) / (1.987*temp))

      conc = self.hybCond['conc']
      target = self.hybCond['target']

      # At equilibrium ([A]t-[AB])([B]t-[AB]) = k[AB].
      # If [A]t >> [B]t, then [AB]/[B]t = [A]t / (k + [A]t).
      if target == 0:
         return conc / (conc + k)
      # Otherwise, solve the (second degree) equation.
      else:
         delta = (k+conc+target)**2 - 4*conc*target
         return (k+conc+target - sqrt(delta)) / (2*target)




##############################################
##                                          ##
##      CLASSES FOR VALIDATOR OBJECT        ##
##                                          ##
##############################################


# Exception declaration for validation.
class ValidationException(Exception):
   """
   Super class for failed validation.
   """
   pass

class oligoTooShortException(ValidationException):
   pass

class oligoTooLongException(ValidationException):
   pass

class oligoSelfDimerException(ValidationException):
   pass

class TmTooLowException(ValidationException):
   pass

class TmTooHighException(ValidationException):
   pass

class DeltaTmTooHighException(ValidationException):
   pass

class HeteroAnnealException(ValidationException):
   pass


# Validator classes.
class Validator(object):
   """
   Abstract (not implemented) validator class.
   """

   DEFAULT_RULES = {}

   def __init__(self, rules={}):
      """Instantiate with a set of rules for validating."""

      # Copy default rules and overwite with constructor argument.
      self.rules = self.DEFAULT_RULES.copy()
      self.rules.update(rules)

   def __call__(self, *args, **kwargs):
      raise NotImplemented


class ConsistencyValidator(Validator):
   """
   Validate that an oligo solution is consistent with a set
   of rules for oligo length, Tm and homodimer formation.
   """

   # See module globals.
   DEFAULT_RULES = STD_RULES

   def __call__(self, oligoSol):
      # Validate length.
      length = len(oligoSol.oligo.seq)
      if length < self.rules['min_length']:
         raise oligoTooShortException(oligoSol)
      if length > self.rules['max_length']:
         raise oligoTooLongException(oligoSol)
      # Validate dimer formation.
      if len(oligoSol.oligo.misprime()) > self.rules['max_self_anneal']:
         raise oligoSelfDimerException(oligoSol)
      # Validate Tm.
      if oligoSol.Tm < self.rules['min_Tm']:
         raise TmTooLowException(oligoSol)
      if oligoSol.Tm > self.rules['max_Tm']:
         raise TmTooHighException(oligoSol)
      if oligoSol.oligo.lanneal() > self.rules['max_lanneal']:
         raise oligoSelfDimerException(oligoSol)
      return


class CompatibilityValidator(Validator):
   """
   Validate that several oligos in solution are compatible with
   a set of rules for Tm difference and heterodimer formation.
   """

   # See module globals.
   DEFAULT_RULES = STD_RULES

   def __call__(self, oligoSolSet):
      # Validate Tm difference.
      TmList = [oligoSol.Tm for oligoSol in oligoSolSet]
      if (max(TmList) - min(TmList)) > self.rules['max_delta_Tm']:
         raise DeltaTmTooHighException(oligoSolSet)
      # Valide hetero dimers.
      oSS = set(oligoSolSet).copy()
      oligo_1 = oSS.pop().oligo
      while oSS:
         for oligo_2 in [oS.oligo for oS in oSS]:
            # Compute the lengths of hetero annealing.
            m1 = len(oligo_1.misprime(oligo_2))
            m2 = len(oligo_2.misprime(oligo_1))
            if max(m1,m2) > self.rules['max_hetero_anneal']:
               raise HeteroAnnealException(oligoSolSet)
            m3 = oligo_1.lanneal(oligo_2)
            if m3 > self.rules['max_lanneal']:
               raise HeteroAnnealException(oligoSolSet)
         oligo_1 = oSS.pop().oligo
      return




##############################################
##                                          ##
##        CLASSES FOR FINDER OBJECTS        ##
##                                          ##
##############################################

class ConsistentOligoFinder(object):
   """
   Find oligos matching a set of consistency rules.
   """

   def __init__(self, hybCond={}, validator=ConsistencyValidator()):
      self.validator = validator
      # Hybridization conditions passed to OligoSol.
      self.hybCond = STD_HYBCOND.copy()
      self.hybCond.update(hybCond)

   def __call__(self, template, extraLeft='', extraRight=''):
      """
      Find consistent oligos that amplify a given template.
      Does not check compatibility.

      RETURN:
        A list of validated 'OligoSol' instances that amplify the
        specified PCR template.
      ARGUMENTS:
        'template': a string containing the sequence to be amplified.
        'extraNtPair': pair of extra nucleotides to be added to
          the 'Oligo' instances.
      """

      # Here function (called twice).
      def get_all_forward_consistent(template, extraNt=''):
         consistent_OligoSols = []
         for i in range(10, len(template)+1):
            sol = OligoSol(Oligo(template[:i], extraNt), self.hybCond)
            try:
               self.validator(sol)
            except oligoTooLongException, TmTooHighException:
               break
            except ValidationException:
               continue
            else:
               consistent_OligoSols.append(sol)
         return consistent_OligoSols

      length = self.validator.rules.get('max_length', 60)
      left_template = template[:length]
      right_template = DNAseq(template[-length:]).revcomp

      return (
         get_all_forward_consistent(left_template, extraLeft),
         get_all_forward_consistent(right_template, extraRight)
      )


class CompatibleOligoFinder(object):
   """
   Find oligos matching a set of compatibility rules.
   """

   def __init__(self, hybCond={}, validator=CompatibilityValidator()):
      self.validator = validator
      self.hybCond = STD_HYBCOND.copy()
      self.hybCond.update(hybCond)

   def oligo_Tm_walk(self, oligoSet):
      """Yield a n-tuple of 'OligoSol' instances from a set of
      n lists of 'OligoSol' instances.

      The oligo with the lowest Tm determines which of the n lists
      is updated.

      This generator allows to perform an ascending Tm search for
      compatible oligo sets (more than two for multiplex PCR)."""

      # Initialize current n-index.
      n = len(oligoSet)
      index = [0] * n
      try:
         current_list = [oligoSet[i][index[i]] for i in range(n)]
         yield current_list
         while True:
            Tm_list = [oS.Tm for oS in current_list]
            (min_Tm,argmin) = min(zip(Tm_list, range(n)))
            index[argmin] += 1
            current_list = [oligoSet[i][index[i]] for i in range(n)]
            yield current_list
      # Stop iteration when the oligo with lowest Tm cannot
      # be updated.
      except IndexError:
         return

   def __call__(self, consistent_OligoSols):
      compatible_sets = []
      # Generator 'oligo_set_walk' simplifies the search.
      for oligo_set in self.oligo_Tm_walk(consistent_OligoSols):
         try:
            self.validator(oligo_set)
         except ValidationException:
            continue
         else:
            compatible_sets.append(oligo_set)

      return compatible_sets



##############################################
##                                          ##
##                FUNCTIONS                 ##
##                                          ##
##############################################


##############################################
###########     primer_search      ###########
##############################################

def score(pair):
   """
   Scoring function for primer pairs. The lower the score
   the better the primer pair.

   RETURN:
      A numeric score.
   ARGUMENT:
      'pair': a tuple '(a,b)', where 'a' and 'b' are instances
            of class 'OligoSol'.
   """

   L,R = pair

   # Add the absolute difference of the Tm to 60C and the
   # absolute difference in Tm between the oligos.
   S1 = abs(L.Tm-333.15) + abs(R.Tm-333.15) + abs(L.Tm - R.Tm)
   # Add the lengths of the self and hetero annealing.
   S2 = len(L.oligo.misprime() + R.oligo.misprime() +
         L.oligo.misprime(R) + R.oligo.misprime(L))
   # Add the lengths of the unspecific annealing.
   S3 = L.oligo.lanneal() + R.oligo.lanneal() + 2 * R.oligo.lanneal(L)

   return S1 + S2 + (.6 * S3)

def chew(n, d=-1):
   """
   A simple iterator to progressively delete the ends of a
   sequence. This is useful to search primer pairs with
   approximate positions.

   RETURN:
      Iterator with start and end position of a subsequence.
      The lengths of the subsequence are nonincreasing.
   ARGUMENT:
      'n': the length of the sequence.
      'd': the maximum length to delete.
   """

   if d < 0 or d > n: d = n-1
   for delta in range(0,d+1):
      for s in range(0, delta+1):
         yield (s,n-delta+s)
   return


def primer_search(template, extraL='', extraR='', **kwargs):
   """
   Perform a product-centric primer search. Return a list of
   consistent and compatible primers to amplify the template.
   Does not perform multiplex search.

   RETURN:
     List of pair of consistent and compatible instances of
     'OligoSol' representing oligos that can amplified the
     given template.
   ARGUMENTS:
     'template': the PCR template as a string.
     'extraL': nucleotides to add to left primer.
     'extraR': nucleotides to add to right primer.
     'kwargs': additional arguments for hybridization
       conditions and validation rules.
   """

   ## PROCESS PARAMETERS ##
   rules = STD_RULES.copy()
   hybCond = STD_HYBCOND.copy()
   rules.update(**kwargs)
   hybCond.update(**kwargs)

   # Instantiate finders.
   f = ConsistentOligoFinder(hybCond, ConsistencyValidator(rules))
   g = CompatibleOligoFinder(hybCond, CompatibilityValidator(rules))

   return g(f(template, extraL, extraR))


##############################################
###########          main          ###########
##############################################


class HitList(list):
   forced = False


def parse_sequences(inputfile):
   """
   Read sequences from input file and return a dictionary.
   """

   seq = dict()
   for line in inputfile:
      if line[0] == '>':
         header = line.rstrip()
         seq[header] = ''
         continue
      seq[header] += line.rstrip()

   return seq


def std_search(seq, args, **kwargs):

   # Update with additional keyword arguments.
   args.__dict__.update(kwargs)

   # Standard search: chew the ends of the sequence until a
   # pair is found. Do not chew at all if 'approx' is equal to 0.

   found = HitList()
   for (s,e) in chew(len(seq), args.approx):
      try:
         primers = primer_search(seq[s:e], **args.__dict__)
         found.extend(primers)
      # Caught non DNA characters: skip.
      except ForbiddenCharacter: pass

   return found



def qPCR_search(seq, args, **kwargs):

   # Using (more restrictive) rules for qPCR.
   rules = qPCR_RULES.copy()
   rules.update(args.__dict__)
   args.__dict__.update(rules)

   # Update with additional keyword arguments.
   args.__dict__.update(kwargs)

   # qPCR search: try all targets between 110 and 180 bp.

   found = HitList()
   for tlen in range(110, 181):
      for s in range(len(seq)-tlen+1):
         e = s+tlen
         try:
            primers = primer_search(seq[s:e], args.__dict__)
            found.extend(primers)
         # Caught non DNA characters: skip.
         except ForbiddenCharacter: pass

   return found



def batch_search(args):
   """
   Wrapper for 'primer_search()'.
   """

   # Input sequences are stored in a dictionary. The key is
   # the title, and the value is the sequence.
   seq = parse_sequences(args.infile) 

   # Set the search function.
   search = qPCR_search if args.qPCR else std_search


   # Valid pairs are stored in a dictionary. The key is the
   # title of the sequence, and the value is a list of
   # pairs of valid 'OligoSol' objects.
   pairs = defaultdict(HitList)
   for title,sequence in seq.items():
      found = search(sequence, args)
      if not found and args.force:
         # When force finding primers, use 'RELAXED_RULES'.
         # TODO: specify that the first search was not successful.
         pairs[title] = search(sequence, args, **RELAXED_RULES)
         pairs[title].forced = True
      else:
         pairs[title] = found

   return pairs


def parse_arguments():
   """
   Parse command line arguments using the 'argparse' module.
   """

   import argparse

   parser = argparse.ArgumentParser()

   # Input-output option.
   parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                   default=sys.stdin)

   # Oligo design options.
   parser.add_argument('--extraL', metavar='L', dest='extraL', default='',
                   type=str, help='extra nucleotides on left primer')

   parser.add_argument('--extraR', metavar='R', dest='extraR', default='',
                   type=str, help='extra nucleotides on right primer')


   # Search options.
   parser.add_argument('--qPCR',  dest='qPCR', action='store_true',
                   help='search primers for qPCR')

   parser.add_argument('--force',  dest='force', action='store_true',
                   help='force find a primer pair')

   parser.add_argument('--approx', metavar='d', dest='approx', default=0,
                   type=int, help='allows d nucleotides shift')


   # Hybridization conditions.
   parser.add_argument('--conc', metavar='C', dest='conc',
         type=float, default=argparse.SUPPRESS,
         help='oligo concentration (%1.f uM)' % 
            (1000000 * STD_HYBCOND['conc']))

   parser.add_argument('--mono', metavar='C', dest='mono',
         type=float, default=argparse.SUPPRESS,
         help='monovalent cations concentration (%.1f mM)' % \
               (1000 * STD_HYBCOND['mono']))

   parser.add_argument('--di', metavar='C', dest='di',
         type=float, default=argparse.SUPPRESS,
         help='concentration of divalent cations (%.1f mM)' % \
               (1000 * STD_HYBCOND['di']))

   # Consistency rules
   parser.add_argument('--min_length', metavar='L', dest='min_length',
         type=int, default=argparse.SUPPRESS,
         help='minimum oligo length (%d)' % STD_RULES['min_length'])

   parser.add_argument('--max_length', metavar='L', dest='max_length',
         type=int, default=argparse.SUPPRESS,
         help='maximum oligo length (%d)' % \
               STD_RULES['max_length'])

   parser.add_argument('--min_Tm', metavar='T', dest='min_Tm',
         type=int, default=argparse.SUPPRESS,
         help='minimum oligo solution Tm (%d)' % \
               (STD_RULES['min_Tm'] + ABS0))

   parser.add_argument('--max_Tm', metavar='T', dest='max_Tm',
         type=int, default=argparse.SUPPRESS,
         help='maximum oligo solution Tm (%d)' % \
               (STD_RULES['max_Tm'] + ABS0))

   parser.add_argument('--max_self_anneal', metavar='M',
         dest='max_self_anneal', default=argparse.SUPPRESS,
         type=int, help='maximum self annealing (%d)' % \
               STD_RULES['max_self_anneal'])

   parser.add_argument('--max_lanneal', metavar='M',
         dest='max_lanneal', default=argparse.SUPPRESS,
         type=int, help='longest unspecific annealing (%d)' % \
               STD_RULES['max_lanneal'])

   # Compatibility rules
   parser.add_argument('--max_delta_Tm', metavar='T',
         dest='max_delta_Tm', type=int, default=argparse.SUPPRESS,
         help='maximum difference between Tm (%d)' % \
               STD_RULES['max_delta_Tm'])

   parser.add_argument('--max_het_anneal', metavar='T',
         dest='max_hetero_anneal', default=argparse.SUPPRESS,
         type=int, help='maximum hetero annealing (%d)' % \
               STD_RULES['max_hetero_anneal'])

   args = parser.parse_args()

   # Post-process arguments.
   if hasattr(args, 'conc'):        args.conc /= 1000000
   if hasattr(args, 'mono'):        args.mono /= 1000
   if hasattr(args, 'di'):          args.di /= 1000
   if hasattr(args, 'max_Tm'):      args.max_Tm -= ABS0
   if hasattr(args, 'min_Tm'):      args.min_Tm -= ABS0

   return args



def show_primers(dict_of_HitList):
   """
   Function to display the results.
   """

   # Just a convenient macro for readability.
   write = sys.stdout.write

   for title,hits in dict_of_HitList.items():
      # Ignore sequences for which no primer pair was found.
      if not hits: continue

      # Rank pairs and show the best (lowest score).
      left,right = min(hits, key=score)

      display_title = title + ' (forced)' if hits.forced else title
      write('\n' + display_title + '\n\n')
      write('(%.1f deg) %s\n' % (left.Tm-273.15, left.oligo))
      write('(%.1f deg) %s\n' % (right.Tm-273.15, right.oligo))
      write('---\n')



if __name__ == '__main__':
   # Parse arguments.
   args = parse_arguments()

   # Search primer pairs.
   dict_of_HitList = batch_search(args)
   
   try:
      show_primers(dict_of_HitList)
   except IOError as e:
      # Piping the output to `head` breaks the pipe. No big deal,
      # we suppress the error so that the user does not panic.
      if e.errno == errno.EPIPE: pass

