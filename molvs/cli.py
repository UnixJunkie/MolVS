# -*- coding: utf-8 -*-
"""
molvs.cli
~~~~~~~~~

This module contains a command line interface for standardization.

:copyright: Copyright 2016 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import argparse
import logging
from rdkit import Chem
import sys

from molvs import Standardizer, Validator


log = logging.getLogger(__name__)


FILETYPES = ['smi', 'mol', 'sdf']


def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            words = line.split('\t')
            smile = words[0]
            name = words[1]
            yield Chem.MolFromSmiles(smile), name


# 'ascii' because this source code file is in utf-8 ?!
name_prop = '_Name'.encode('ascii')


def MySDMolSupplier(filename):
    for mol in Chem.SDMolSupplier(filename):
        name = mol.GetProp(name_prop)
        yield mol, name


def reader_of_file(fn):
    if fn.endswith('.smi'):
        return RobustSmilesMolSupplier(fn)
    elif fn.endswith('.sdf'):
        return MySDMolSupplier(fn)
    else:
        return None

def write_smiles_out(out, mol, name):
    out.write(Chem.MolToSmiles(mol))
    if name != None:
        if name.endswith('\n'):
            out.write('\t%s' % name)
        else:
            out.write('\t%s\n' % name)


def write_sdf_out(out, mol, name):
    mol.SetProp(name_prop, name)
    out.write(mol)


def writer_of_file(fn):
    if fn.endswith('.smi'):
        return open(fn, 'w')
    elif fn.endswith('.sdf'):
        return Chem.SDWriter(fn)
    else:
        return None


class MolvsParser(argparse.ArgumentParser):

    def error(self, message):
        sys.stderr.write('Error: %s\n\n'.encode() % message)
        self.print_help()
        sys.exit(2)


def main():
    """Main function for molvs command line interface."""

    # Root options
    parser = MolvsParser(epilog='use "molvs <command> -h" to show help for a specific command')
    subparsers = parser.add_subparsers(title='Available commands')

    # Options common to all commands

    common_parser = MolvsParser(add_help=False)
    common_parser.add_argument('infile', nargs='?', help='input filename', type=argparse.FileType('r'), default=sys.stdin)
    common_parser.add_argument('-i', '--intype', help='input filetype', choices=FILETYPES)
    common_parser.add_argument('-:', '--smiles', help='input SMILES instead of file', metavar='<smiles>')
    common_parser.add_argument('-O', '--outfile', help='output filename', type=argparse.FileType('w'), default=sys.stdout, metavar='<outfile>')

    # Standardize options
    standardize_parser = subparsers.add_parser('standardize', help='standardize a molecule', parents=[common_parser])
    standardize_parser.add_argument('-o', '--outtype', help='output filetype', choices=FILETYPES)
    standardize_parser.set_defaults(func=standardize_main)

    # Validate options
    validate_parser = subparsers.add_parser('validate', help='validate a molecule', parents=[common_parser])
    validate_parser.set_defaults(func=validate_main)

    args = parser.parse_args()
    try:
        args.func(args)
    except Exception as e:
        sys.stderr.write('Error: %s\n\n'.encode() % e.message)
        parser.print_help()
        sys.exit(2)


def _read_mol(args):
    if args.smiles:
        return Chem.MolFromSmiles(args.smiles)
    elif args.intype in {'smi', 'smiles'} or args.infile.name.endswith('smi') or args.infile.name.endswith('smiles'):
        return Chem.MolFromSmiles(args.infile.read())
    elif args.intype in {'mol', 'sdf'} or args.infile.name.endswith('mol') or args.infile.name.endswith('sdf'):
        return Chem.MolFromMolBlock(args.infile.read())
    else:
        return Chem.MolFromSmiles(args.infile.read())


def _write_mol(mol, args):
    if args.outtype in {'smi', 'smiles'} or args.outfile.name.endswith('smi') or args.outfile.name.endswith('smiles'):
        args.outfile.write(Chem.MolToSmiles(mol))
        args.outfile.write('\n')
    elif args.outtype in {'mol', 'sdf'} or args.outfile.name.endswith('mol') or args.outfile.name.endswith('sdf'):
        args.outfile.write(Chem.MolToMolBlock(mol))
    else:
        args.outfile.write(Chem.MolToSmiles(mol))
        args.outfile.write('\n')


def standardize_main(args):
    in_fn = args.infile.name
    out_fn = args.outfile.name
    reader = reader_of_file(in_fn)
    writer = writer_of_file(out_fn)
    s = Standardizer()
    if (reader is not None) and (writer is not None):
        smiles_out = out_fn.endswith('.smi')
        sdf_out = out_fn.endswith('.sdf')
        for mol, name in reader:
            std = s.standardize(mol)
            if smiles_out:
                write_smiles_out(writer, std, name)
            elif sdf_out:
                write_sdf_out(writer, std, name)
            else:
                assert(False)
        writer.close()
    else:
        mol = _read_mol(args)
        mol = s.standardize(mol)
        _write_mol(mol, args)


def validate_main(args):
    mol = _read_mol(args)
    v = Validator()
    logs = v.validate(mol)
    for log in logs:
        args.outfile.write(log)
        args.outfile.write('\n')
