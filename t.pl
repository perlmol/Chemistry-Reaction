#!/usr/bin/perl

use strict;
use warnings;
use Chemistry::File::SMILES;
use Chemistry::Pattern;
use lib 'lib';
use Chemistry::Reaction;

open F, "reactions.txt" or die;
my @reacts;
while (<F>) {
    chomp;
    my ($subst, $prod, $name) = split " ", $_, 3;
    my $s = Chemistry::Pattern->parse($subst, format=>'smiles');
    my $p  = Chemistry::Pattern->parse($prod, format=>'smiles');
    my %m;
    for (my $i = 1; $i <= $s->atoms; $i++) {
        $m{$s->atoms($i)} = $p->atoms($i);
    }
    my $react = Chemistry::Reaction->new($s, $p, \%m, name => $name);
    push @reacts, $react;
}

my @mols = Chemistry::Mol->read('mols.smi');
for my $react (@reacts) {
    for my $mol (@mols) {
        my @products;
        printf "%s\t%s\t", $react->name, $mol->name;
        my $subst = $react->substrate;
        while ($subst->match($mol)) {
            my $new_mol = $mol->clone;
            my @map = map { $new_mol->by_id($_) } $subst->atom_map;
            $react->forward($new_mol, @map);
            push @products, $new_mol;
        }
        print join(", ", map {$_->print(format=>'smiles')} @products), "\n";
    }
}
