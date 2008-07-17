#
# = bio/data/aa.rb - Amino Acids
#
# Copyright::	Copyright (C) 2001, 2005
#		Toshiaki Katayama <k@bioruby.org>
# License::	The Ruby License
#
# $Id: aa.rb,v 0.22 2007/04/06 04:44:51 k Exp $
#

module Bio

class AminoAcid

  module Data

    # IUPAC code
    # * http://www.iupac.org/
    # * http://www.chem.qmw.ac.uk/iubmb/newsletter/1999/item3.html
    # * http://www.ebi.ac.uk/RESID/faq.html

    NAMES = {

      'A' => 'Ala',
      'C' => 'Cys',
      'D' => 'Asp',
      'E' => 'Glu',
      'F' => 'Phe',
      'G' => 'Gly',
      'H' => 'His',
      'I' => 'Ile',
      'K' => 'Lys',
      'L' => 'Leu',
      'M' => 'Met',
      'N' => 'Asn',
      'P' => 'Pro',
      'Q' => 'Gln',
      'R' => 'Arg',
      'S' => 'Ser',
      'T' => 'Thr',
      'V' => 'Val',
      'W' => 'Trp',
      'Y' => 'Tyr',
      'B' => 'Asx',	# D/N
      'Z' => 'Glx',	# E/Q
      'J' => 'Xle',	# I/L
      'U' => 'Sec',	# 'uga' (stop)
      'O' => 'Pyl',	# 'uag' (stop)
      'X' => 'Xaa',	# (unknown)
     
      'Ala' => 'alanine',
      'Cys' => 'cysteine',
      'Asp' => 'aspartic acid',
      'Glu' => 'glutamic acid',
      'Phe' => 'phenylalanine',
      'Gly' => 'glycine',
      'His' => 'histidine',
      'Ile' => 'isoleucine',
      'Lys' => 'lysine',
      'Leu' => 'leucine',
      'Met' => 'methionine',
      'Asn' => 'asparagine',
      'Pro' => 'proline',
      'Gln' => 'glutamine',
      'Arg' => 'arginine',
      'Ser' => 'serine',
      'Thr' => 'threonine',
      'Val' => 'valine',
      'Trp' => 'tryptophan',
      'Tyr' => 'tyrosine',
      'Asx' => 'asparagine/aspartic acid [DN]',
      'Glx' => 'glutamine/glutamic acid [EQ]',
      'Xle' => 'isoleucine/leucine [IL]',
      'Sec' => 'selenocysteine',
      'Pyl' => 'pyrrolysine',
      'Xaa' => 'unknown [A-Z]',

    }
    
    ## Taken from http://matrixscience.com/help/aa_help.html
    # The data in this table are for amino acid residues. To 
    # calculate the mass of a neutral peptide or protein, sum 
    # the residue masses plus the masses of the terminating groups
    # (e.g. H at the N-terminus and OH at the C-terminus).
    MONOISOTOPICMASS = {
      'A' => 71.03712 ,
      'R' => 156.10112,
      'N' => 114.04293,
      'D' => 115.02695,
      'C' => 103.00919,
      'E' => 129.0426 ,
      'Q' => 128.05858,
      'G' => 57.02147 ,
      'H' => 137.05891,
      'I' => 113.08407,
      'L' => 113.08407,
      'K' => 128.09497,
      'M' => 131.04049,
      'F' => 147.06842,
      'P' => 97.05277 ,
      'S' => 87.03203 ,
      'T' => 101.04768,
      'U' => 150.95364,
      'W' => 186.07932,
      'Y' => 163.06333,
      'V' => 99.06842
    }
    AVGMASS = {
      'A' => 71.08 ,
      'R' => 156.19,
      'N' => 114.1 ,
      'D' => 115.09,
      'C' => 103.14,
      'E' => 129.12,
      'Q' => 128.13,
      'G' => 57.05 ,
      'H' => 137.14,
      'I' => 113.16,
      'L' => 113.16,
      'K' => 128.17,
      'M' => 131.19,
      'F' => 147.18,
      'P' => 97.12 ,
      'S' => 87.08 ,
      'T' => 101.1 ,
      'U' => 150.03,
      'W' => 186.21,
      'Y' => 163.18,
      'V' => 99.13
    }

    def weight(x = nil, monoisotopic=true)
      table = monoisotopic ? MONOISOTOPICMASS : AVGMASS
      unless x 
        return table
      end
      total = 0.0
      x.each_byte do |byte|
        aa = byte.chr.upcase
        if table[aa]
          total += table[aa]
        else
          raise "Error: invalid amino acid '#{aa}'"
        end
      end
      total += NucleicAcid.weight[:water]
    end


    def [](x)
      NAMES[x]
    end

    # backward compatibility
    def names
      NAMES
    end
    alias aa names

    def name(x)
      str = NAMES[x]
      if str and str.length == 3
        NAMES[str]
      else
        str
      end
    end

    def to_1(x)
      case x.to_s.length
      when 1
        x
      when 3
        three2one(x)
      else
        name2one(x)
      end
    end
    alias one to_1

    def to_3(x)
      case x.to_s.length
      when 1
        one2three(x)
      when 3
        x
      else
        name2three(x)
      end
    end
    alias three to_3

    def one2three(x)
      if x and x.length != 1
        raise ArgumentError
      else
        NAMES[x]
      end
    end

    def three2one(x)
      if x and x.length != 3
        raise ArgumentError
      else
        reverse[x]
      end
    end

    def one2name(x)
      if x and x.length != 1
        raise ArgumentError
      else
        three2name(NAMES[x])
      end
    end

    def name2one(x)
      str = reverse[x.to_s.downcase]
      if str and str.length == 3
        three2one(str)
      else
        str
      end
    end

    def three2name(x)
      if x and x.length != 3
        raise ArgumentError
      else
        NAMES[x]
      end
    end

    def name2three(x)
      reverse[x.downcase]
    end

    def to_re(seq)
      replace = {
        'B' => '[DNB]',
        'Z' => '[EQZ]',
        'J' => '[ILJ]',
        'X' => '[ACDEFGHIKLMNPQRSTVWYUOX]',
      }
      replace.default = '.'

      str = seq.to_s.upcase
      str.gsub!(/[^ACDEFGHIKLMNPQRSTVWYUO]/) { |aa|
        replace[aa]
      }
      Regexp.new(str)
    end


    private


    def reverse
      hash = Hash.new
      NAMES.each do |k, v|
        hash[v] = k
      end
      hash
    end

  end


  # as instance methods
  include Data

  # as class methods
  extend Data


  private


  # override when used as an instance method to improve performance
  alias orig_reverse reverse
  def reverse
    unless @reverse
      @reverse = orig_reverse
    end
    @reverse
  end

end

end # module Bio


if __FILE__ == $0

  puts "### aa = Bio::AminoAcid.new"
  aa = Bio::AminoAcid.new

  puts "# Bio::AminoAcid['A']"
  p Bio::AminoAcid['A']
  puts "# aa['A']"
  p aa['A']

  puts "# Bio::AminoAcid.name('A'), Bio::AminoAcid.name('Ala')"
  p Bio::AminoAcid.name('A'), Bio::AminoAcid.name('Ala')
  puts "# aa.name('A'), aa.name('Ala')"
  p aa.name('A'), aa.name('Ala')

  puts "# Bio::AminoAcid.to_1('alanine'), Bio::AminoAcid.one('alanine')"
  p Bio::AminoAcid.to_1('alanine'), Bio::AminoAcid.one('alanine')
  puts "# aa.to_1('alanine'), aa.one('alanine')"
  p aa.to_1('alanine'), aa.one('alanine')
  puts "# Bio::AminoAcid.to_1('Ala'), Bio::AminoAcid.one('Ala')"
  p Bio::AminoAcid.to_1('Ala'), Bio::AminoAcid.one('Ala')
  puts "# aa.to_1('Ala'), aa.one('Ala')"
  p aa.to_1('Ala'), aa.one('Ala')
  puts "# Bio::AminoAcid.to_1('A'), Bio::AminoAcid.one('A')"
  p Bio::AminoAcid.to_1('A'), Bio::AminoAcid.one('A')
  puts "# aa.to_1('A'), aa.one('A')"
  p aa.to_1('A'), aa.one('A')

  puts "# Bio::AminoAcid.to_3('alanine'), Bio::AminoAcid.three('alanine')"
  p Bio::AminoAcid.to_3('alanine'), Bio::AminoAcid.three('alanine')
  puts "# aa.to_3('alanine'), aa.three('alanine')"
  p aa.to_3('alanine'), aa.three('alanine')
  puts "# Bio::AminoAcid.to_3('Ala'), Bio::AminoAcid.three('Ala')"
  p Bio::AminoAcid.to_3('Ala'), Bio::AminoAcid.three('Ala')
  puts "# aa.to_3('Ala'), aa.three('Ala')"
  p aa.to_3('Ala'), aa.three('Ala')
  puts "# Bio::AminoAcid.to_3('A'), Bio::AminoAcid.three('A')"
  p Bio::AminoAcid.to_3('A'), Bio::AminoAcid.three('A')
  puts "# aa.to_3('A'), aa.three('A')"
  p aa.to_3('A'), aa.three('A')

  puts "# Bio::AminoAcid.one2three('A')"
  p Bio::AminoAcid.one2three('A')
  puts "# aa.one2three('A')"
  p aa.one2three('A')

  puts "# Bio::AminoAcid.three2one('Ala')"
  p Bio::AminoAcid.three2one('Ala')
  puts "# aa.three2one('Ala')"
  p aa.three2one('Ala')

  puts "# Bio::AminoAcid.one2name('A')"
  p Bio::AminoAcid.one2name('A')
  puts "# aa.one2name('A')"
  p aa.one2name('A')

  puts "# Bio::AminoAcid.name2one('alanine')"
  p Bio::AminoAcid.name2one('alanine')
  puts "# aa.name2one('alanine')"
  p aa.name2one('alanine')

  puts "# Bio::AminoAcid.three2name('Ala')"
  p Bio::AminoAcid.three2name('Ala')
  puts "# aa.three2name('Ala')"
  p aa.three2name('Ala')

  puts "# Bio::AminoAcid.name2three('alanine')"
  p Bio::AminoAcid.name2three('alanine')
  puts "# aa.name2three('alanine')"
  p aa.name2three('alanine')

  puts "# Bio::AminoAcid.to_re('BZACDEFGHIKLMNPQRSTVWYU')"
  p Bio::AminoAcid.to_re('BZACDEFGHIKLMNPQRSTVWYU')

end

