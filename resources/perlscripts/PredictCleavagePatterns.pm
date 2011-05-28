package PredictCleavagePatterns;

our @EXPORT = ("getTetramerCleavageValues", "getHeptamerCleavageValues", "predictSequenceCleavagePattern", "predictSequenceCleavagePatternWithTetramers", "predictSequenceCleavagePatternWithTrimers", "getTrimerCleavageValues", "getPentamerCleavageValues", "getMinTetramerScore", "getMinTrimerScore", "getScrambledHeptamerCleavageValues");
use base "Exporter";
use Exporter;

use strict;
use warnings;

my $defaultTetramerCleavageValues = getTetramerCleavageValues();
my $defaultHeptamerCleavageValues = getHeptamerCleavageValues($defaultTetramerCleavageValues);
my $minTetramerScore = getMinTetramerScore($defaultTetramerCleavageValues);


my $defaultTrimerCleavageValues = getTrimerCleavageValues();
my $defaultPentamerCleavageValues = getPentamerCleavageValues($defaultTrimerCleavageValues);
my $minTrimerScore = getMinTrimerScore($defaultTrimerCleavageValues);

#####################
#### SUBROUTINES ####
#####################


sub predictSequenceCleavagePattern
{
	my $sequence = uc shift(@_);
	return predictSequenceCleavagePatternWithTetramers($sequence);
}

sub getMinTetramerScore
{
	if(scalar(@_) == 0 && defined($minTetramerScore))
	{
		return $minTetramerScore;
	}
	
	my $tetramers = scalar(@_) > 0 ? shift(@_) : $defaultTetramerCleavageValues;

	my $min = 99999999999;
	foreach my $tetramerScores (values %$tetramers)
	{
		foreach my $tetramerScore (@$tetramerScores)
		{
			$min = $tetramerScore < $min ? $tetramerScore : $min;
		}
	}
	return $min;
}

sub predictSequenceCleavagePatternWithTetramers
{
	my $sequence = uc shift(@_);
	my $tetramerCleavageValues = scalar(@_) > 0 ? shift(@_) : $defaultTetramerCleavageValues;
	my $heptamerCleavageValues = scalar(@_) > 0 ? shift(@_) : $defaultHeptamerCleavageValues;
	
	my $seqLength = length($sequence);
	my $cleavagePattern = [];

	
	# iterate over sequence, computing cleavage pattern at each nucleotide, saving it into an array
	foreach my $seqIndex (0..length($sequence)-1)
	{
		

		my $heptamerStart = $seqIndex-3;
		my $heptamerEnd = $seqIndex +3;
		if($heptamerStart >= 0 && $heptamerEnd < length($sequence) )
		{
			push(@$cleavagePattern, $heptamerCleavageValues->{substr($sequence,$heptamerStart,7)});
		}
		else
		{
			# Because prediction is based on tetramer we want to average over all tetramers that contain nucleotide.
			# However, first two nucleotides and last two do not have three tetramers like the rest.
			# So, save all windows that contain nucleotide and save its position in the window.
			my $tetramerWindows = {};
			foreach my $windowStart ($seqIndex-3..$seqIndex)
			{
				if($windowStart >= 0 and $windowStart+3 < length($sequence))
				{
					my $positionInWindow = $seqIndex-$windowStart;
					$tetramerWindows->{$positionInWindow} = substr($sequence,$windowStart,4);
				}
			}
		
			# average the predicted values for each window
			my $meanCleavageValue = 0;
			foreach my $positionInWindow (keys %$tetramerWindows)
			{
				$meanCleavageValue = $meanCleavageValue + $tetramerCleavageValues->{$tetramerWindows->{$positionInWindow}}->[$positionInWindow];
			}
			$meanCleavageValue = $meanCleavageValue/scalar(keys %$tetramerWindows);

			push(@$cleavagePattern, $meanCleavageValue);
		}
	}
	return $cleavagePattern;
}


sub generateAllNmers
{
	my $seqLength = shift(@_);
	my $sequences = scalar(@_) > 0 ? shift(@_) : [ "" ];

	my $newSequences = [];
	foreach my $sequence (@$sequences)
	{
		push(@$newSequences, "A$sequence");
		push(@$newSequences, "T$sequence");
		push(@$newSequences, "C$sequence");
		push(@$newSequences, "G$sequence");
	}
	
	if($seqLength > 1)
	{
		$sequences = generateAllNmers($seqLength-1, $newSequences);
	}
	else
	{
		$sequences = $newSequences;	
	}
	return $sequences;
}


sub getHeptamerCleavageValues
{	
	my $tetramerCleavageValues = scalar(@_) > 0 ? shift(@_) : getTetramerCleavageValues();

	my $heptamerCleavageValues = {};
	my $heptamers = generateAllNmers(7);

	foreach my $heptamer (@$heptamers)
	{
		my $tetramerWindows = {3=>substr($heptamer,0,4), 2=>substr($heptamer,1,4), 1=>substr($heptamer,2,4), 0=>substr($heptamer,3,4) };
		my $meanCleavageValue = 0;
		foreach my $positionInWindow (keys %$tetramerWindows)
		{
			$meanCleavageValue = $meanCleavageValue + $tetramerCleavageValues->{$tetramerWindows->{$positionInWindow}}->[$positionInWindow];
		}
		$meanCleavageValue = $meanCleavageValue/scalar(keys %$tetramerWindows);
		
		$heptamerCleavageValues->{$heptamer} = $meanCleavageValue;
	}
	return $heptamerCleavageValues;	
}


sub getScrambledHeptamerCleavageValues
{
	my $scrambledNucs = scalar(@_) > 0 ? shift(@_) : {"A"=>"C", "T"=>"G", "C"=>"T", "G"=>"A"};
	

	my $scrambledHeptamerValues = {};
	my $heptamerCleavageValues = getHeptamerCleavageValues();
	foreach my $heptamerSeq (keys %$heptamerCleavageValues)
	{
		my @splitHeptamer = split("", $heptamerSeq);
		my @newHeptamer = ();
		foreach my $oldNuc (@splitHeptamer)
		{
			push(@newHeptamer, $scrambledNucs->{$oldNuc});
		}
		my $newHeptamerSeq = join("", @newHeptamer);
		$scrambledHeptamerValues->{$newHeptamerSeq} = $heptamerCleavageValues->{$heptamerSeq};
	}
	
	return $scrambledHeptamerValues;	

}


sub getTetramerCleavageValues
{
	#colums = 
	#	tetramer, 
	#	instances in database,
	#	mean at position 1, 
	#	mean at position 2, 
	#	mean at position 3,
	#	mean at position 4, 
	#	st. dev at position 1, 
	#	st. dev at position 2, 
	#	st. dev at position 3,
	#	st. dev at position 4 
	my @cleavageData  = qw(
		GAGA         8    0.742    0.058    1.117    0.025  0.240  0.391  0.246  0.288
		TACA         8    0.440   -0.185    1.078   -0.473  0.645  0.351  0.377  0.139
		CTGG         8    0.417   -0.227    0.427    1.309  0.656  0.222  0.218  0.355
		CTAC         8   -0.062   -0.377   -0.119    0.684  0.930  0.564  0.275  0.429
		AGAT         8   -0.196    1.129   -0.001    0.489  0.313  0.286  0.338  0.178
		GATA         8    0.730    0.151    0.801   -0.547  0.526  0.528  0.338  0.282
		TCAG         8    0.565    0.392   -0.432    1.370  0.646  0.465  0.285  0.502
		ATTT        10   -0.176    0.396   -0.544   -0.977  0.601  0.818  0.564  0.427
		TAGC         8    0.176   -0.577    1.025    0.791  0.798  0.255  0.391  0.648
		GACT         8    0.837    0.143    0.536   -0.371  0.314  0.368  0.490  0.157
		CCAG         8    0.264    0.169   -0.407    1.777  0.485  0.471  0.438  0.348
		ATCC         8    0.341    0.780    0.042   -0.147  0.506  0.623  0.348  0.267
		TCGC         8    0.235   -0.120    0.712    1.229  0.574  0.204  0.147  0.378
		GTCC         8    0.964    0.752    0.399    0.057  0.993  0.514  0.273  0.359
		TACC         8    0.856   -0.115    0.814    0.093  1.292  0.358  0.428  0.178
		GTCG         8    0.786    0.611    0.044    0.583  0.348  0.708  0.366  0.641
		GCGA         8    1.053    1.233    0.823    0.802  0.387  0.428  0.146  0.286
		TATT         8    0.599   -0.366    0.564   -0.110  0.883  0.265  0.321  0.181
		ATCA         8   -0.302    0.475    0.122   -0.491  0.699  1.349  0.796  0.452
		AGTA         8   -0.167    0.924    0.601   -0.349  0.417  0.554  0.746  0.222
		GACA         8    0.258    0.116    0.979   -0.421  0.235  0.415  0.376  0.618
		CAGA         8    0.700   -0.645    1.557    0.263  0.464  0.269  0.185  0.231
		GGGT         8    0.163    0.792    0.930    0.179  0.439  0.280  0.270  0.436
		TCGT         8    0.444    0.175    0.703    1.713  0.496  0.558  0.360  0.252
		CGGT         8    0.574    0.594    1.912    1.119  0.888  0.351  0.566  0.431
		CGTT         8    0.348    1.090    2.109    0.676  0.734  0.345  0.784  0.362
		CGTG         8    0.516    0.696    1.746    0.821  0.773  0.220  0.461  0.463
		AGGT         8   -0.185    1.115    1.282    0.355  0.344  0.690  0.608  0.301
		CTTC         8    0.684    0.328    0.447    0.365  0.570  0.344  0.245  0.362
		TCCT         8    0.675    0.185   -0.060    0.269  0.504  0.351  0.211  0.474
		TGTG         8    0.656    0.724    0.941    0.871  0.780  0.193  0.217  0.495
		TAAG         8    0.859   -0.216    0.106    0.717  0.715  0.223  0.232  0.308
		TCGG         8   -0.028   -0.132    0.135    1.318  0.321  0.460  0.354  0.430
		CCTG         8    0.456    0.027   -0.061    0.745  0.575  0.180  0.231  0.240
		TTCG         8    0.129   -0.073   -0.273    0.490  0.494  0.358  0.296  0.203
		AACT         8   -0.013    0.158   -0.093   -0.934  0.390  0.516  0.445  0.278
		GTTG        10    1.156    1.354   -0.116    0.410  0.616  0.954  0.414  0.336
		CCCT        10    0.387   -0.194    0.088    0.309  0.609  0.302  0.308  0.517
		TCGA         8    0.594    0.126    0.779    0.637  0.605  0.571  0.420  0.462
		ACCG         8    0.048    0.705    0.008    0.582  0.571  0.279  0.337  0.218
		TGAA         8    0.162    0.563    0.060   -0.045  0.876  0.602  0.457  0.200
		CGAA         8    0.336    0.542    0.543    0.250  0.346  0.468  0.345  0.235
		ATAC         8   -0.212    0.916   -0.520    0.516  0.762  0.761  0.246  0.319
		TCCG         8    0.363    0.257    0.050    0.846  0.979  0.266  0.274  0.395
		CTGT         8    0.651    0.005    1.041    1.263  0.953  0.293  0.252  0.397
		TAAA         8    0.092   -0.358   -0.068   -0.318  0.973  0.709  0.475  0.325
		GGCT         8    0.557    1.023    1.143    0.088  0.167  0.248  0.307  0.123
		GTCA         8    1.375    1.015    0.687   -0.061  0.611  0.582  0.230  0.208
		GTAA         8    0.681    0.675   -0.313    0.222  0.629  0.733  0.159  0.249
		TAGA        10    0.550   -0.413    0.980    0.109  0.557  0.434  0.351  0.296
		ATGC         8   -0.152    0.721    0.316    0.809  0.714  1.130  0.509  0.316
		GCTG         8    0.912    1.295   -0.144    0.861  0.444  0.438  0.307  0.442
		TAAC         8    0.179   -0.393    0.204    0.561  1.078  0.333  0.306  0.343
		GGTG         8    0.616    1.576    0.739    0.444  0.395  0.561  0.517  0.676
		TTGG         8    0.354   -0.187    0.064    1.082  0.468  0.409  0.316  0.295
		TGGT         8    0.264    0.162    1.410    0.602  0.621  0.184  0.288  0.396
		GGGA         8    0.676    0.727    0.800   -0.024  0.613  0.232  0.331  0.385
		GGCG         8    0.153    1.085    1.021    0.476  0.417  0.542  0.485  0.556
		TATC         8    0.504   -0.348    0.843    0.158  0.599  0.319  0.611  0.459
		GGCC         8    0.567    0.840    0.634   -0.352  0.170  0.715  0.379  0.337
		ACGA         8   -0.226    0.633    0.521    0.440  0.336  0.499  0.501  0.373
		ACGG         8   -0.383    0.811    0.518    1.086  0.494  0.325  0.259  0.304
		GAAA         8    0.649    0.204    0.046   -0.685  0.487  0.329  0.363  0.265
		TGTT         8    0.180    0.829    1.267    0.016  0.562  0.325  0.506  0.328
		TACT         8    0.258   -0.359    0.523    0.042  0.825  0.339  0.531  0.232
		TCTT         8    0.072   -0.268    0.052    0.462  0.604  0.565  0.429  0.266
		TTGA         8   -0.517   -0.773    0.365   -0.035  0.621  0.585  0.313  0.209
		GCCT         8    0.810    0.992   -0.014    0.118  0.550  0.407  0.363  0.663
		TTAC         8   -0.320   -0.328   -0.416    0.644  0.750  0.553  0.267  0.469
		GGAA         8    0.566    1.044    0.155   -0.081  0.266  0.257  0.323  0.335
		TCTC         8    0.370   -0.005   -0.108    0.426  0.598  0.409  0.335  0.310
		ATAT         8   -0.197    0.882   -0.515    0.306  0.539  0.531  0.399  0.267
		CGCC         8    0.613    0.412    0.977   -0.036  0.676  0.280  0.560  0.366
		ATTC         8    0.291    0.836   -0.054   -0.441  0.623  0.603  0.624  0.230
		GTAG         8    0.871    1.005   -0.287    1.119  0.592  0.617  0.191  0.177
		CAAA         8    0.061   -0.418    0.359   -0.018  0.995  0.564  0.411  0.160
		ACAC         8    0.083    1.133   -0.362    1.315  0.510  0.480  0.258  0.164
		GTAT         8    0.716    0.944   -0.262    0.941  0.295  0.634  0.427  0.386
		TATA         8    0.813   -0.060    0.974   -0.461  0.576  0.683  0.769  0.348
		GATC         8    1.031    0.298    0.510   -0.111  0.340  0.445  0.216  0.343
		ATCT         8   -0.016    0.633   -0.287   -0.416  0.485  0.332  0.393  0.332
		CGAG         8    0.630    0.731    0.494    1.177  0.676  0.396  0.388  0.284
		TAGT        10    0.496   -0.360    1.113    0.785  0.326  0.439  0.348  0.205
		GGGC        10    0.344    0.778    0.741    0.798  0.298  0.370  0.484  0.437
		TCCC         8   -0.043   -0.104   -0.332    0.089  0.416  0.288  0.273  0.311
		CTCA         8    0.140   -0.287    0.485   -0.214  0.413  0.555  0.328  0.347
		GTAC         8    0.771    0.960   -0.169    1.140  0.629  0.867  0.462  0.429
		AAGG         8   -0.385   -0.239    0.275    0.389  0.287  0.329  0.473  0.342
		CGGA         8    0.161    0.262    1.028    0.252  0.418  0.351  0.191  0.314
		ACAA         8   -0.520    0.673   -0.666    0.192  0.663  1.046  0.423  0.390
		GAGG         8    0.575    0.132    0.757    1.090  0.372  0.301  0.266  0.486
		GGTC         8    0.970    1.429    0.417    0.146  0.783  0.663  0.658  0.519
		GTTT         8    0.987    1.166   -0.092   -0.289  0.390  0.571  0.414  0.282
		CCCG         8    0.434   -0.037    0.083    0.896  0.824  0.338  0.292  0.464
		CCAA         8   -0.293   -0.114   -0.141    0.693  0.607  0.446  0.184  0.331
		GAGC         8    0.936   -0.005    0.902    0.400  0.536  0.342  0.439  0.326
		TCTG         8    0.213    0.119   -0.136    1.023  0.327  0.471  0.438  0.462
		CCCA         8   -0.238   -0.330    0.482    0.125  0.404  0.581  0.274  0.213
		ACGT         8   -0.195    0.603    0.731    1.349  0.453  0.311  0.612  0.494
		AGCA         8   -0.170    1.123    0.635   -0.572  0.342  0.500  0.670  0.209
		CTTA         8    0.275   -0.110    0.456   -0.192  0.423  0.672  0.512  0.250
		GGGG         8    0.056    0.722    0.646    0.656  0.743  0.463  0.364  0.240
		AAGC         8    0.140    0.246    0.619    0.288  0.388  0.226  0.203  0.459
		AGCC         8    0.127    1.337    0.416    0.010  0.406  0.553  0.518  0.359
		CCGA         8    0.632    0.078    0.863    0.621  0.287  0.169  0.464  0.401
		TTTC         8    0.611   -0.252   -0.366   -0.582  0.860  0.585  0.368  0.390
		CCGT         8    0.510    0.089    1.117    2.249  0.667  0.319  0.317  0.303
		CATA         8    0.223   -0.703    1.251   -0.218  0.386  0.240  0.419  0.163
		CGCT         8    0.192    0.354    1.498    0.031  0.505  0.163  0.293  0.295
		GTTA         8    0.950    0.701    0.067   -0.413  0.274  0.971  0.561  0.223
		AGTG         8   -0.458    0.903    0.525    0.501  0.346  0.708  0.530  0.353
		CTCG         8    0.437    0.179    0.506    0.889  0.391  0.310  0.373  0.282
		AATT         8   -0.088   -0.018   -0.044   -0.867  0.498  0.236  0.409  0.205
		CGAC         8    0.518    0.694    0.582    0.972  0.917  0.409  0.421  0.282
		AGCG         8   -0.243    0.949    0.747    0.492  0.390  0.514  0.539  0.447
		TTCA         8    0.825    0.295    0.172   -0.372  1.061  0.758  0.667  0.331
		TGCA         8    0.384    0.499    1.396   -0.129  1.227  0.722  0.725  0.804
		CCCC         8    0.271   -0.268   -0.117    0.076  0.666  0.349  0.331  0.411
		GGAG        10    0.468    0.797   -0.106    0.803  0.418  0.297  0.315  0.308
		GTGA         8    0.932    0.835    0.896    0.508  0.875  0.751  0.416  0.188
		AGGC         8   -0.518    0.600    0.699    0.802  0.512  0.306  0.333  0.359
		TCTA        10    0.574    0.103   -0.124   -0.225  0.656  0.672  0.590  0.368
		ATTA         8   -0.131    0.344   -0.235   -0.696  0.274  0.604  0.545  0.303
		CAGG         8    0.512   -0.764    1.329    0.987  1.082  0.358  0.548  0.457
		AGAA        10   -0.382    1.096    0.020   -0.200  0.447  0.474  0.548  0.458
		AATC         8    0.009    0.025   -0.170   -0.693  0.530  0.900  0.897  0.370
		GTGC         8    1.363    1.291    1.053    1.528  0.540  0.616  0.395  0.267
		CTCT         8    0.311   -0.062    0.475    0.330  0.700  0.364  0.370  0.273
		GCAA         8    0.827    1.265   -0.138    0.342  0.521  0.483  0.747  0.293
		AGAG         8   -0.120    0.893   -0.012    0.906  0.512  0.396  0.306  0.405
		AATA         8   -0.569   -0.519    0.010   -0.960  0.480  0.740  0.756  0.376
		ACTA        10   -0.278    0.263   -0.410   -0.264  0.276  0.774  0.611  0.329
		ATGT         8   -0.022    1.006    0.600    0.764  0.523  0.694  0.337  0.325
		AGGG        10   -0.383    0.655    0.690    0.597  0.586  0.477  0.260  0.457
		AGAC         8   -0.278    1.094   -0.016    0.639  0.461  0.523  0.163  0.334
		TGGG         8    0.143   -0.148    0.748    0.921  0.305  0.352  0.346  0.328
		GGCA         8    0.700    1.070    1.013   -0.534  0.536  0.400  0.233  0.208
		CATC         8    0.958   -0.075    1.232    0.295  0.608  0.110  0.409  0.324
		CCTA         8    0.169    0.001    0.594    0.281  0.308  0.361  0.566  0.289
		TCAA         8    0.691    0.734   -0.002    0.621  0.730  0.240  0.278  0.175
		CACC         8    0.623   -0.286    0.946   -0.020  0.527  0.404  0.413  0.618
		CTGC         8    0.884    0.030    1.126    1.585  0.551  0.400  0.221  0.229
		GCTT         8    0.762    1.077    0.086    0.421  0.345  0.452  0.293  0.263
		CTAT         8    0.458    0.154    0.020    1.346  0.701  0.578  0.445  0.489
		CTTT         8    0.125    0.009    0.421    0.171  1.001  0.385  0.266  0.199
		TAGG         8    0.211   -0.805    0.573    0.922  1.019  0.483  0.271  0.365
		ATCG         8   -0.123    0.527   -0.228    0.368  0.367  0.473  0.425  0.246
		CGGC         8    0.691    0.229    1.414    1.128  0.605  0.440  0.236  0.195
		GTCT         8    1.089    0.866    0.285    0.222  0.365  0.610  0.355  0.393
		AACA         8   -0.243   -0.106    0.371   -0.904  0.547  0.876  0.864  0.289
		CACT         8    0.771   -0.177    1.079    0.083  0.555  0.307  0.374  0.562
		GCAC         8    0.436    0.669   -0.679    0.998  0.661  0.634  0.221  0.320
		TGAG         8    0.214    0.563    0.206    0.950  0.861  0.380  0.318  0.288
		AACG         8   -0.267    0.134    0.506    0.114  0.690  0.292  0.209  0.248
		GCCC        10    0.818    0.534   -0.193    0.118  0.514  0.666  0.409  0.471
		TCAT         8    0.709    0.502   -0.240    1.468  0.635  0.354  0.336  0.522
		GCGG         8    0.573    0.867   -0.117    1.493  0.363  0.714  0.402  0.559
		TAAT         8    0.205   -0.358   -0.030    0.399  0.446  0.398  0.353  0.342
		AAAT        10   -0.068   -0.127   -0.767   -0.884  0.392  0.594  0.661  0.451
		AGGA         8   -0.709    0.400    0.613   -0.250  0.562  0.561  0.468  0.266
		TGCT         8    0.782    1.002    1.453    0.225  0.404  0.312  0.250  0.225
		TGGA         8    0.204    0.174    0.970    0.117  0.996  0.396  0.339  0.295
		CGCG        12    0.315    0.601    1.447    0.487  0.374  0.474  0.284  0.540
		GCCG         8    0.978    0.742   -0.170    0.913  0.864  0.432  0.380  0.536
		GAAC         8    0.588    0.244    0.070    0.302  0.714  0.456  0.146  0.316
		AAAA         8   -0.215   -0.095   -0.536   -1.126  0.650  0.352  0.443  0.490
		CCTC         8    0.386   -0.202    0.107    0.597  0.809  0.449  0.443  0.447
		ATAG         8   -0.158    0.713   -0.588    0.899  0.332  0.402  0.333  0.317
		CAGT         8    0.686   -0.330    1.859    1.245  0.413  0.213  0.289  0.238
		TTTA         8    0.073   -0.342   -0.426   -0.849  0.666  0.659  0.593  0.421
		GGAT         8    0.294    0.996    0.038    0.386  0.408  0.294  0.450  0.259
		AAGA         8   -0.058   -0.108    0.639   -0.397  0.426  0.338  0.281  0.305
		TTGT         8    0.777    0.070    0.596    0.996  0.671  0.345  0.102  0.351
		CCGC         8    0.326   -0.110    0.672    1.462  0.405  0.285  0.386  0.598
		GAAT         8    1.064    0.264   -0.068    0.098  0.270  0.458  0.382  0.322
		CCAT         8    0.573    0.188   -0.236    1.755  0.691  0.263  0.360  0.454
		TGCG         8    0.204    0.898    1.363    0.627  0.603  0.318  0.477  0.331
		TGTA         8    0.574    0.633    1.065   -0.258  0.208  0.395  0.534  0.137
		CTCC        10    0.437   -0.394    0.207   -0.062  0.942  0.366  0.347  0.311
		GTTC         8    1.314    1.254    0.247    0.071  0.589  0.919  0.598  0.407
		AGCT         8   -0.308    0.712    0.588   -0.215  0.457  0.260  0.232  0.267
		AAAG         8   -0.021    0.070   -0.324    0.162  0.388  0.149  0.376  0.347
		CATG         8    0.768   -0.155    1.373    0.724  0.462  0.368  0.719  0.468
		CACA         8    0.137   -0.511    1.400   -0.134  0.842  0.405  0.362  0.120
		TTAG         8    0.276   -0.137   -0.800    0.627  1.106  0.683  0.489  0.425
		CAAG         8    0.829   -0.411    0.260    0.734  0.673  0.300  0.290  0.313
		AATG        10    0.057   -0.281    0.110   -0.024  0.577  0.569  0.913  0.430
		CGGG         8    0.034    0.037    1.027    1.054  0.631  0.487  0.309  0.303
		TGCC         8    0.598    0.713    1.182   -0.118  1.036  0.319  0.175  0.160
		GCTA         8    0.635    1.162    0.134    0.083  0.348  0.456  0.265  0.174
		TCAC         8   -0.466   -0.161   -0.465    1.009  1.259  0.768  0.414  0.418
		CTTG         8    0.299   -0.171    0.108    0.620  0.695  0.384  0.387  0.132
		TTAA        12    0.305    0.056   -0.514   -0.013  0.713  0.710  0.384  0.233
		TGAT         8   -0.196    0.826    0.443    0.774  0.904  0.548  0.392  0.259
		GACG        10    0.468    0.131    0.622    0.459  0.309  0.205  0.314  0.340
		CGTC        10    0.841    0.781    1.411    0.536  0.442  0.731  0.708  0.289
		ACTT         8   -0.007    0.558   -0.175    0.103  0.492  0.547  0.655  0.326
		GCCA         8    0.664    0.773   -0.096   -0.407  0.298  0.654  0.301  0.350
		CACG         8    0.383   -0.560    1.152    0.552  0.915  0.319  0.222  0.426
		TTTG         8    0.507   -0.105   -0.404    0.315  1.004  0.844  0.702  0.377
		GAGT         8    0.675    0.329    1.072    0.690  0.230  0.548  0.332  0.568
		GCAG         8    0.938    1.318   -0.477    1.721  0.566  0.637  0.240  0.428
		TTTT         8    0.111   -0.212   -0.648   -0.809  0.626  0.523  0.551  0.392
		GTGT         8    0.747    0.957    0.650    1.232  0.380  0.507  0.149  0.338
		GATT         8    1.138    0.557    0.520   -0.381  0.217  0.530  0.327  0.282
		TTCC         8    0.967    0.222   -0.058   -0.024  0.728  0.485  0.265  0.284
		CGCA         8    0.320    0.542    1.351   -0.392  0.622  0.330  0.431  0.189
		AAAC         8   -0.379   -0.158   -0.668   -0.445  0.696  0.795  0.613  0.499
		TTGC         8    0.504   -0.073    0.617    1.471  0.764  0.331  0.166  0.467
		ACAT         8   -0.096    1.117   -0.224    0.937  0.487  0.456  0.494  0.268
		CGAT         8    0.586    1.019    0.881    0.829  0.554  0.258  0.344  0.272
		ACCA         8   -0.305    0.457   -0.008   -0.310  0.406  0.667  0.492  0.355
		ACTC         8    0.302    0.423   -0.508    0.071  0.393  0.371  0.439  0.195
		ATAA         8   -0.563    0.526   -0.562   -0.314  0.886  1.137  0.600  0.422
		AACC        10   -0.155   -0.088    0.077   -0.576  0.546  0.630  0.581  0.291
		ACTG         8   -0.226    0.652   -0.200    0.687  0.299  0.723  0.370  0.391
		ATGG         8   -0.378    0.559    0.067    0.815  0.442  0.410  0.382  0.429
		GGAC         8    0.222    0.574   -0.001    0.519  0.551  0.420  0.327  0.397
		TGTC         8    0.627    0.701    0.983    0.391  0.912  0.214  0.166  0.289
		GAAG         8    1.134    0.272    0.051    0.535  0.316  0.432  0.243  0.286
		CAAT        10    0.714   -0.113    0.425    0.659  0.572  0.657  0.350  0.396
		CCGG         8    0.200   -0.086    0.585    1.485  0.336  0.461  0.364  0.630
		TGGC         8    0.402    0.412    1.210    1.168  0.573  0.297  0.400  0.399
		AGTT         8    0.081    1.358    0.801   -0.313  0.348  0.508  0.786  0.382
		CATT        10    0.916   -0.218    1.329    0.308  0.511  0.317  0.475  0.392
		ATGA         8   -0.263    0.177    0.527   -0.006  0.436  0.591  0.614  0.341
		ACCT         8    0.072    0.681   -0.216   -0.059  0.447  0.442  0.481  0.331
		CCAC         8    0.608    0.273   -0.027    1.255  0.675  0.356  0.311  0.468
		CGTA         8    0.526    0.415    1.553    0.101  0.543  0.353  0.449  0.325
		TGAC         8   -0.288    0.558    0.207    0.683  0.664  0.385  0.138  0.502
		ACGC         8   -0.040    0.820    0.337    1.208  0.460  0.340  0.292  0.487
		GACC         8    1.266    0.382    0.679    0.125  0.354  0.426  0.302  0.266
		ATTG        10   -0.343    0.326   -0.529    0.181  0.459  0.605  0.583  0.453
		GATG         8    1.071    0.355    0.647    0.362  0.193  0.564  0.260  0.384
		GTGG         8    0.856    0.868    0.041    1.133  0.171  0.657  0.416  0.471
		GCAT         8    1.034    1.143   -0.331    1.121  0.398  0.549  0.249  0.327
		GCTC         8    0.781    1.149    0.052    0.591  0.329  0.515  0.197  0.221
		CTAA         8    0.125   -0.013    0.008    0.273  0.587  0.839  0.303  0.142
		TTAT         8    0.162    0.094   -0.531    0.417  0.288  0.321  0.280  0.272
		AAGT         8    0.017    0.196    0.615   -0.194  0.275  0.295  0.434  0.285
		TTCT        10    0.163   -0.215   -0.557   -0.383  0.745  0.479  0.590  0.220
		CTGA         8    0.141   -0.348    0.721    0.449  0.385  0.279  0.431  0.294
		TACG         8   -0.382   -0.566    0.570    0.835  0.432  0.209  0.337  0.535
		CTAG         8    0.557    0.037   -0.229    1.220  0.541  0.547  0.467  0.379
		GGTA         8    0.510    1.068    0.365   -0.524  0.263  0.610  0.397  0.242
		TCCA         8    0.321    0.250    0.138   -0.220  0.698  0.426  0.291  0.250
		GCGT         8    1.129    1.431    0.740    1.809  0.382  0.217  0.260  0.688
		CCTT         8    0.668    0.016    0.093    0.445  0.514  0.379  0.529  0.541
		GGTT        10    0.518    1.365    0.653   -0.208  0.480  0.372  0.441  0.181
		AGTC         8    0.365    1.453    0.625    0.327  0.530  0.314  0.504  0.280
		ACCC         8    0.051    0.625   -0.144    0.235  0.668  0.818  0.559  0.214
		CAGC         8    0.885   -0.257    1.576    0.908  0.474  0.335  0.435  0.286
		GCGC         8    0.624    0.900    0.442    1.522  0.616  0.274  0.364  0.231
		ACAG         8   -0.152    0.904   -0.679    1.452  0.548  0.319  0.424  0.274
		TATG         8    0.158   -0.514    0.629    0.610  0.331  0.292  0.374  0.277
		CAAC        10    0.862   -0.169    0.678    0.598  0.515  0.338  0.232  0.453	
	);

	my $fieldsPerLine = 10;
	my $tetramerCleavageValues = { };
	while(scalar(@cleavageData) > 0)
	{
		my @nextLineFields = ();
		foreach(1..$fieldsPerLine)
		{
			push(@nextLineFields, shift(@cleavageData));
		}
		$tetramerCleavageValues->{$nextLineFields[0]} = [ $nextLineFields[2], $nextLineFields[3], $nextLineFields[4], $nextLineFields[5] ];
	}
	return $tetramerCleavageValues;
}









sub getMinTrimerScore
{
	if(scalar(@_) == 0 && defined($minTrimerScore))
	{
		return $minTrimerScore;
	}
	
	my $trimers = scalar(@_) > 0 ? shift(@_) : $defaultTrimerCleavageValues;

	my $min = 99999999999;
	foreach my $trimerScores (values %$trimers)
	{
		foreach my $trimerScore (@$trimerScores)
		{
			$min = $trimerScore < $min ? $trimerScore : $min;
		}
	}
	return $min;
}

sub predictSequenceCleavagePatternWithTrimers
{
	my $sequence = uc shift(@_);
	my $trimerCleavageValues = scalar(@_) > 0 ? shift(@_) : $defaultTrimerCleavageValues;
	my $pentamerCleavageValues = scalar(@_) > 0 ? shift(@_) : $defaultPentamerCleavageValues;
	
	my $seqLength = length($sequence);
	my $cleavagePattern = [];

	
	# iterate over sequence, computing cleavage pattern at each nucleotide, saving it into an array
	foreach my $seqIndex (0..length($sequence)-1)
	{
		my $pentamerStart = $seqIndex-2;
		my $pentamerEnd = $seqIndex +2;
		if($pentamerStart >= 0 && $pentamerEnd < length($sequence) )
		{
			push(@$cleavagePattern, $pentamerCleavageValues->{substr($sequence,$pentamerStart,5)});
		}
		else
		{
			# Because prediction is based on trimer we want to average over all trimers that contain nucleotide.
			# However, first two nucleotides and last two do not have three trimers like the rest.
			# So, save all windows that contain nucleotide and save its position in the window.
			my $trimerWindows = {};
			foreach my $windowStart ($seqIndex-2..$seqIndex)
			{
				if($windowStart >= 0 and $windowStart+2 < length($sequence))
				{
					my $positionInWindow = $seqIndex-$windowStart;
					$trimerWindows->{$positionInWindow} = substr($sequence,$windowStart,3);
				}
			}
		
			# average the predicted values for each window
			my $meanCleavageValue = 0;
			foreach my $positionInWindow (keys %$trimerWindows)
			{
				$meanCleavageValue = $meanCleavageValue + $trimerCleavageValues->{$trimerWindows->{$positionInWindow}}->[$positionInWindow];
			}
			$meanCleavageValue = $meanCleavageValue/scalar(keys %$trimerWindows);

			push(@$cleavagePattern, $meanCleavageValue);
		}
	}
	return $cleavagePattern;
}



sub getPentamerCleavageValues
{
	my $trimerCleavageValues = scalar(@_) > 0 ? shift(@_) : getTrimerCleavageValues();

	my $pentamerCleavageValues = {};
	my $pentamers = generateAllNmers(5);

	foreach my $pentamer (@$pentamers)
	{
		#print "$pentamer\t";
		my $trimerWindows = {2=>substr($pentamer,0,3), 1=>substr($pentamer,1,3), 0=>substr($pentamer,2,3) };
		my $meanCleavageValue = 0;
		foreach my $positionInWindow (keys %$trimerWindows)
		{
			$meanCleavageValue = $meanCleavageValue + $trimerCleavageValues->{$trimerWindows->{$positionInWindow}}->[$positionInWindow];
		}
		$meanCleavageValue = $meanCleavageValue/scalar(keys %$trimerWindows);
		#print "$meanCleavageValue\n";
		
		$pentamerCleavageValues->{$pentamer} = $meanCleavageValue;
	}
	return $pentamerCleavageValues;
}



sub getTrimerCleavageValues
{

	#colums = 
	#	trimer, 
	#	instances in database,
	#	mean at position 1, 
	#	mean at position 2, 
	#	mean at position 3, 
	#	st. dev at position 1, 
	#	st. dev at position 2, 
	#	st. dev at position 3 
	my @cleavageData  = qw(
		AAA 34 -0.165 -0.081 -0.585 0.535 0.519 0.547
		AAC 36 -0.147 0.064 0.227 0.529 0.615 0.593
		AAG 32 -0.072 0.024 0.537 0.387 0.352 0.379
		AAT 38 -0.104 -0.159 -0.003 0.556 0.645 0.735
		ACA 32 -0.171 0.957 -0.483 0.573 0.636 0.437
		ACC 34 -0.030 0.597 -0.123 0.518 0.557 0.466
		ACG 34 -0.195 0.707 0.488 0.429 0.368 0.456
		ACT 34 -0.066 0.462 -0.328 0.421 0.625 0.531
		AGA 36 -0.284 1.056 -0.001 0.441 0.411 0.356
		AGC 32 -0.148 1.030 0.597 0.417 0.505 0.504
		AGG 34 -0.445 0.690 0.813 0.526 0.562 0.487
		AGT 34 -0.060 1.162 0.641 0.492 0.554 0.611
		ATA 32 -0.283 0.759 -0.546 0.653 0.740 0.396
		ATC 32 -0.025 0.604 -0.088 0.556 0.767 0.526
		ATG 34 -0.156 0.656 0.392 0.549 0.772 0.485
		ATT 38 -0.111 0.517 -0.336 0.527 0.697 0.588
		CAA 36 0.636 -0.263 0.444 0.736 0.491 0.348
		CAC 32 0.479 -0.383 1.144 0.736 0.379 0.373
		CAG 32 0.696 -0.499 1.580 0.648 0.356 0.417
		CAT 34 0.728 -0.284 1.298 0.559 0.360 0.497
		CCA 32 0.288 0.129 -0.203 0.694 0.402 0.349
		CCC 36 0.215 -0.213 0.115 0.658 0.393 0.356
		CCG 32 0.417 -0.007 0.809 0.458 0.323 0.422
		CCT 34 0.380 -0.042 0.168 0.591 0.354 0.509
		CGA 32 0.517 0.747 0.625 0.634 0.410 0.389
		CGC 36 0.355 0.491 1.333 0.535 0.352 0.427
		CGG 32 0.365 0.280 1.345 0.684 0.441 0.501
		CGT 34 0.574 0.748 1.687 0.626 0.514 0.653
		CTA 36 0.360 0.013 -0.055 0.732 0.633 0.369
		CTC 34 0.338 -0.156 0.406 0.655 0.448 0.362
		CTG 32 0.523 -0.135 0.829 0.695 0.332 0.396
		CTT 32 0.346 0.014 0.358 0.702 0.483 0.380
		GAA 34 0.828 0.184 -0.030 0.518 0.466 0.360
		GAC 34 0.693 0.189 0.699 0.481 0.356 0.393
		GAG 34 0.749 0.131 0.950 0.368 0.394 0.339
		GAT 32 0.992 0.340 0.619 0.365 0.514 0.300
		GCA 32 0.809 1.099 -0.406 0.567 0.609 0.452
		GCC 34 0.818 0.747 -0.123 0.570 0.561 0.359
		GCG 36 0.852 1.178 0.517 0.489 0.506 0.468
		GCT 32 0.773 1.171 0.032 0.365 0.450 0.277
		GGA 34 0.392 0.849 0.014 0.425 0.357 0.352
		GGC 34 0.521 1.000 0.963 0.403 0.477 0.391
		GGG 36 0.340 0.770 0.796 0.558 0.330 0.374
		GGT 34 0.646 1.359 0.550 0.525 0.556 0.509
		GTA 32 0.760 0.896 -0.258 0.532 0.696 0.326
		GTC 34 1.073 0.892 0.361 0.628 0.666 0.368
		GTG 32 0.975 0.988 0.660 0.578 0.634 0.521
		GTT 36 1.117 1.115 -0.004 0.484 0.850 0.492
		TAA 36 0.283 -0.364 0.035 0.867 0.437 0.345
		TAC 32 0.293 -0.306 0.746 0.934 0.352 0.462
		TAG 36 0.377 -0.522 0.937 0.686 0.433 0.387
		TAT 32 0.518 -0.322 0.753 0.644 0.438 0.549
		TCA 32 0.375 0.367 -0.285 0.954 0.579 0.367
		TCC 34 0.297 0.151 -0.045 0.691 0.344 0.302
		TCG 32 0.311 0.012 0.582 0.541 0.470 0.415
		TCT 36 0.255 -0.069 -0.100 0.626 0.595 0.446
		TGA 32 -0.027 0.627 0.229 0.821 0.479 0.358
		TGC 32 0.492 0.778 1.348 0.866 0.473 0.449
		TGG 32 0.253 0.150 1.085 0.645 0.363 0.415
		TGT 32 0.509 0.722 1.064 0.667 0.288 0.394
		TTA 36 0.128 -0.064 -0.559 0.771 0.601 0.378
		TTC 34 0.500 0.041 -0.201 0.837 0.555 0.549
		TTG 36 0.405 -0.245 0.372 0.946 0.567 0.373
		TTT 34 0.343 -0.235 -0.505 0.780 0.618 0.565
	);

		
	my $fieldsPerLine = 8;
	my $trimerCleavageValues = { };
	while(scalar(@cleavageData) > 0)
	{
		my @nextLineFields = ();
		foreach(1..$fieldsPerLine)
		{
			push(@nextLineFields, shift(@cleavageData));
		}
		$trimerCleavageValues->{$nextLineFields[0]} = [ $nextLineFields[2], $nextLineFields[3], $nextLineFields[4] ];
	}
	return $trimerCleavageValues;
}
