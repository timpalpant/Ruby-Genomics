##
# This script attempts to setup all script dependencies automatically on OS X
# You must first have Xcode installed (available from the App Store on 10.7)
#
# What it installs:
#   - Homebrew (package manager): http://mxcl.github.com/homebrew/
#     - Git (SCM): http://www.fftw.org/
#     - FFTW3 library: http://www.fftw.org/
#   - RVM (Ruby version manager): http://beginrescueend.com/rvm/install/
#     - Ruby 1.9.3: http://www.ruby-lang.org/en/
#   - Bundler gem (Ruby gem manager): http://gembundler.com/
#   - The RubyGenomics scripts: https://github.com/timpalpant/Ruby-Genomics
#     - All required gems
##

# Check to ensure that Xcode is installed
hash xcodebuild 2>&- || { echo >&2 "You must first install Xcode. Aborting."; exit 1; }

# Install Homebrew
/usr/bin/ruby -e "$(curl -fsSL https://raw.github.com/gist/323731)"

# Install FFTW3
brew install fftw

# Install Git
brew install git

# Install RVM (Ruby Version Manager)
bash < <(curl -s https://raw.github.com/wayneeseguin/rvm/master/binscripts/rvm-installer )
echo '[[ -s "$HOME/.rvm/scripts/rvm" ]] && . "$HOME/.rvm/scripts/rvm" # Load RVM function' >> ~/.bash_profile
source ~/.bash_profile

# Install Ruby 1.9.3
rvm install 1.9.3
rvm use 1.9.3 --default

# Install Bundler
gem install bundler

# Check out Ruby-Genomics scripts from GitHub
git clone git://github.com/timpalpant/Ruby-Genomics.git

# Install all required gems with Bundler
cd Ruby-Genomics
bundle install

# Run the tests
bundle exec rake unit_test
bundle exec rake test