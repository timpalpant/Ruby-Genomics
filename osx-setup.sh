#!/usr/bin/env bash

##
# This script attempts to setup all script dependencies automatically on OS X
# You must first have Xcode installed (available from the App Store on 10.6.5 and higher)
#
# What it installs:
#   - Homebrew (package manager): http://mxcl.github.com/homebrew/
#     - Git (SCM): http://www.git-scm.com/
#     - FFTW3 library: http://www.fftw.org/
#   - Samtools
#   - Tabix
#   - UCSC Tools
#   - RVM (Ruby version manager): http://beginrescueend.com/rvm/install/
#     - Ruby 1.9.3: http://www.ruby-lang.org/en/
#   - Bundler gem (Ruby gem manager): http://gembundler.com/
#   - The RubyGenomics scripts: https://github.com/timpalpant/Ruby-Genomics
#     - All required gems
##

# Check to ensure that Xcode is installed
hash xcodebuild 2>&- || { echo >&2 "You must first install Xcode. Aborting."; exit 1; }

# Install Homebrew (if it is not already)
hash brew 2>&- || { /usr/bin/ruby -e "$(curl -fsSL https://raw.github.com/gist/323731)"; }

# Install FFTW3
brew install fftw

# Install Git (if it is not already)
hash git 2>&- || { brew install git; }

# Install RVM (Ruby Version Manager)
hash rvm 2>&- || {
  curl -O https://raw.github.com/wayneeseguin/rvm/master/binscripts/rvm-installer
  sh rvm-installer
  rm rvm-installer
  echo '[[ -s "$HOME/.rvm/scripts/rvm" ]] && . "$HOME/.rvm/scripts/rvm" # Load RVM function' >> ~/.bash_profile
  source ~/.bash_profile
}

# Install Ruby 1.9.2
rvm install 1.9.2
rvm use 1.9.2 --default

# Check out Ruby-Genomics scripts from GitHub
git clone git://github.com/timpalpant/Ruby-Genomics.git
cd Ruby-Genomics

# Install Tabix
hash tabix 2>&- || {
  cd ext
  tar xfj tabix-0.2.5.tar.bz2
  cd tabix-0.2.5
  make
  echo 'export PATH=$PATH:'`pwd` >> ~/.bash_profile
  cd ../.. 
}

# Install Samtools
hash samtools 2>&- || {
  cd ext
  tar xfj samtools-0.1.18.tar.bz2
  cd samtools-0.1.18
  make
  echo 'export PATH=$PATH:'`pwd` >> ~/.bash_profile
  cd ../..
}

# Download UCSC binaries
cd ext
mkdir ucsc-binaries
cd ucsc-binaries
echo 'export PATH=$PATH:'`pwd` >> ~/.bash_profile
for b in bedGraphToBigWig bedSort bedToBigBed bigBedInfo bigBedSummary bigBedToBed bigWigInfo bigWigSummary bigWigToBedGraph bigWigToWig faToNib faToTwoBit nibFrag twoBitInfo twoBitToFa wigCorrelate wigToBigWig
do
  hash $b 2>&- || {
    curl -O http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.i386/$b
    chmod a+x $b
  }
done
cd ../..

# Install Bundler
gem install bundler

# Install all required gems with Bundler
gem install fftw3 -- --with-fftw-dir=/usr/local/Cellar/fftw/3.3
bundle install

# Run the tests
source ~/.bash_profile
bundle exec rake unit_test
bundle exec rake test
