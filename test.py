# -*- coding:utf-8 -*-

import oligoz
import unittest

class TestOligoz(unittest.TestCase):
   def test_fasta_search(self):
      print oligoz.fasta_search(open('seq.fasta'))
   def test_fasta_search(self):
      print oligoz.fasta_search(open('seq.fasta'),
            extraL='GATAGG', extraR='GAGAG')

if __name__ == '__main__':
       unittest.main()
