#!/usr/bin/perl

use strict;
use warnings;
use Chemistry::File::SMILES;
use Chemistry::Pattern;
use lib 'lib';
use Chemistry::Reaction;

my $s = Chemistry::Pattern->parse('C=CC=C.C=C', format=>'smiles');
my $p = Chemistry::Pattern->parse('C1C=CCCC1', format=>'smiles');
my %m;
for (my $i = 1; $i <= $s->atoms; $i++) {
    $m{$s->atoms($i)} = $p->atoms($i);
}
my $react = Chemistry::Reaction->new($s, $p, \%m);

my $mol = Chemistry::Mol->parse('C=CC=C.C=C', format=>'smiles');
my $subst = $react->substrate;
while ($subst->match($mol)) {
    my $new_mol = $mol->clone;
    my @map = map { $new_mol->by_id($_) } $subst->atom_map;
    $react->forward($new_mol, @map);
    $new_mol->printf("Product: %s\n");
}

