# picarro.R
# R script (written using 3.0.3) to process Picarro instrument data
# BBL April 2014

# To use:
# 	1. Set SYSTEM_VOLUME and CHAMBER_AREA variables below
#	2. Set MEAS_INTERVAL variable
#	3. If using plot-specific ancillary data, set PLOTDATA variable
#			and make sure file is in correct format
#	4. Set INPUT_DIR to location of analyzer output files
#	5. Check the flux calculations in compute_flux() for your use case
#	6. source() this file.
#			start with a simple test case, and carefully check log file and outputs!

# Important variable definitions, esp. data source & destination
SCRIPTNAME		<- "picarro.R"
#INPUT_DIR		<- "sampledata/"
INPUT_DIR <- "~/Documents/Work/Current/Bailey\ SFA/DWP\ experiment/Data/Data\ from\ Sarah"
OUTPUT_DIR		<- "outputs/"
LOG_DIR			<- "logs/"

# The optional SOLENOIDDATA file *must* have a 'solenoid_valves' field
# It *may* have 'Mass' and/or 'Area' fields, which will be divided into the computed flux
#	- if present, the 'Area' field will override CHAMBER_AREA below
# It *may* have a 'Volume' field, cm3, which will be added to SYSTEM_VOLUME for each plot
SOLENOIDDATA		<- "sampledata/solenoid_valves.csv"


# Picarros output a ton of data. For testing, we may want to subsample
# e.g. 0.1 = sample 10%. Set to 1 for no subsampling
SUBSAMPLE_FRACTION	<- 0.01

SYSTEM_VOLUME	<- 15.15 / 100^3		# m^3
CHAMBER_AREA	<- 5.5					# cm^2

SEPARATOR		<- "-------------------"


# -----------------------------------------------------------------------------
# Time-stamped output function
printlog <- function( msg="", ..., ts=TRUE, cr=TRUE ) {
	if( ts ) cat( date(), " " )
	cat( msg, ... )
	if( cr ) cat( "\n")
} # printlog

# -----------------------------------------------------------------------------
# Print dimensions of data frame
printdims <- function( d, dname=deparse( substitute( d ) ) ) {
	stopifnot( is.data.frame( d ) )
	printlog( dname, "rows =", nrow( d ), "cols =", ncol( d ) )
} # printdims

# -----------------------------------------------------------------------------
# Save a ggplot figure
saveplot <- function( pname, p=last_plot(), ptype=".pdf" ) {
	stopifnot( file.exists( OUTPUT_DIR ) )
	fn <- paste0( OUTPUT_DIR, "/", pname, ptype )
	printlog( "Saving", fn )
	ggsave( fn, p )
} # saveplot

# -----------------------------------------------------------------------------
# Save a data frame
savedata <- function( df, extension=".csv" ) {
	stopifnot( file.exists( OUTPUT_DIR ) )
	fn <- paste0( OUTPUT_DIR, "/", deparse( substitute( df ) ), extension )
	printlog( "Saving", fn )
	write.csv( df, fn, row.names=F )
} # savedata

# -----------------------------------------------------------------------------
# Open a csv file and return data
read_csv <- function( fn, datadir=".", ... ) {
	fqfn <- paste( datadir, fn, sep="/" )
	printlog( "Opening", fqfn )
	stopifnot( file.exists( fqfn ) )
	read.csv( fqfn, stringsAsFactors=F, ... )
} # read_csv

# -----------------------------------------------------------------------------
# Load requested libraries
loadlibs <- function( liblist ) {
	printlog( "Loading libraries..." )
	loadedlibs <- vector()
	for( lib in liblist ) {
		printlog( "Loading", lib )
		loadedlibs[ lib ] <- require( lib, character.only=T )
		if( !loadedlibs[ lib ] )
			warning( "this package is not installed!" )
	}
	invisible( loadedlibs )
} # loadlibs

# -----------------------------------------------------------------------------
# read a process a single output file, returning data frame
read_outputfile <- function( fn ) {
	fqfn <- paste( INPUT_DIR, fn, sep="/" )
	printlog( "Reading", fqfn )
	stopifnot( file.exists( fqfn ) )
	d <- read.table( fqfn, header=T )
	printdims( d )

	if( SUBSAMPLE_FRACTION < 1 ) {
		printlog( "Subsampling at", SUBSAMPLE_FRACTION, "..." )
		d <- d[ sample( nrow( d ), nrow( d ) * SUBSAMPLE_FRACTION ), ]
		printdims( d )
	}
	
	# Add ancillary data
	d$file <- basename( fn )
	d$dir <- dirname( fn )
	
	return( d )
} # read_outputfile


# -----------------------------------------------------------------------------
# read plot data, if it exists, and merge with analyzer data
read_plotdata <- function( fn=PLOTDATA ) {
	d <- NULL
	if( file.exists( fn ) ) {
		d <- read_csv( fn )
		if( any( names( d )=="Plot" ) ) {
			printlog( "Plot data read OK" )
		} else {
			printlog( "Plot data file read, but no 'Plot' field!" )
			warning( "No plot field!" )
		}
	} else {
		printlog( "Plot data file", fn, "not found" )
	}
	
	return( d )
} # read_plotdata

# -----------------------------------------------------------------------------
# compute fluxes
compute_flux <- function( d ) {

	m <- lm( CO2_Ref ~ Sec, data=d )
	resp_raw <- as.numeric( coef( m )[ 2 ] )	# i.e. the slope
	
	# We want to convert raw respiration (d[CO2]/dt) to a flux using
	# A = dC/dt * V/S * Pa/RT (e.g. Steduto et al. 2002), where
	# 	A is CO2 flux (umol/m2/s)
	#	dC/dt is raw respiration as above (mole fraction/s)
	# 	V is total chamber volume (m3)
	#		...correcting for varying headspaces in the cores, if applicable
	#	S is ground surface area (m2), if applicable
	# 	M is sample dry mass (g), if applicable
	#	Pa is atmospheric pressure (kPa)
	#	R is universal gas constant (8.3 x 10-3 m-3 kPa mol-1 K-1)
	#	T is air temperature (K)

	S 			<- CHAMBER_AREA		# note cm2, not m2!
	if( any( names( d )=="Area" ) ) {
		 S <- mean( d$Area )
	}
	V 			<- SYSTEM_VOLUME
	if( any( names( d )=="Volume" ) ) {
		V <- V + mean( d$Volume )
	}
	M 			<- 1.0
	if( any( names( d )=="Mass" ) ) {
		 M <- mean( d$Mass )
	}
	
	Pa 			<- 101						# kPa
	R 			<- 8.3e-3					# m-3 kPa mol-1 K-1
	
	# Calculate mass- (or area-) corrected respiration, umol/g soil/s or umol/cm2/s
	resp_corrected <- resp_raw * V/S/M * Pa/( R*( 273.1+Tair ) )

	# Convert from umol/g soil/s to mgC/kg soil/day or whatever
	# NOTE: you probably want to change this line for your specific setup
	flux <- resp_corrected / 1e6 * 12 * 1000 * 1000 * 60 * 60 * 24
printlog(nrow(d),length(flux))
	return( c( Tair=Tair, V=V, S=S, Day=mean( d$Day ), Month=mean( d$Month ), N=nrow( d ), flux=flux ) )
}

# ==============================================================================
# Main

if( !file.exists( OUTPUT_DIR ) ) {
	printlog( "Creating", OUTPUT_DIR )
	dir.create( OUTPUT_DIR )
}
if( !file.exists( LOG_DIR ) ) {
	printlog( "Creating", LOG_DIR )
	dir.create( LOG_DIR )
}

sink( paste0( LOG_DIR, SCRIPTNAME, ".txt" ), split=T )

printlog( "Welcome to", SCRIPTNAME )

loadlibs( c( "ggplot2", "reshape2", "plyr" ) )
theme_set( theme_bw() )

tf <- tempfile()
printlog( "tempfile is", tf )
alldata <- data.frame()
filelist <- list.files( path=INPUT_DIR, pattern="dat$", recursive=T )
for( f in 1:length( filelist ) ) {
	printlog( SEPARATOR )
	printlog( "Processing file", f, "of", length( filelist ) )
	d <- read_outputfile( filelist[ f ] )
	write.table( d, tf, row.names=F, append=( f>1 ), col.names=( f==1 ), sep="," )
}

printlog( SEPARATOR )
printlog( "All done. Reading data back in..." )
alldata <- read.csv( tf )
printdims( alldata )

# Fractional solenoid values mean that the analyzer was shifting
#  between two samples. Discard these
printlog( "Removing fractional solenoid_valves" )
alldata <- subset( alldata, solenoid_valves==trunc( solenoid_valves ) )
printdims( alldata )

printlog( "** NOTE **" )
printlog( "** Here this script assumes data are stored in a particular way" )
printlog( "** We assume first folder level of path contains treatment info," )
printlog( "** the second is a rep, and there are exactly two levels." )
printlog( "**  i.e. {INPUT_DIR}/treatmentname/repnum/{files}" )
printlog( "** This is very specific to a particular setup-change as necessary." )
printlog( "Splitting file path data..." )
alldata <- cbind( alldata, colsplit( alldata$dir, "/", names=c( "treatment", "rep" ) ) )

if( any( names( alldata )=="solenoid_valves" ) ) {
	sv <- read_csv( "solenoid_valves.csv", INPUT_DIR, comment.char="#" )
	printlog( "Merging Picarro and solenoid_valves data..." )
	alldata <- merge( alldata, sv )
}

	# temporary
	
	printlog( "Getting rid of crazy values..." )
	alldata <- alldata[ alldata$CH4_dry < 5, ]

	printlog( "Focusing on 'real' cores only..." )
	alldata <- alldata[ !is.na( alldata$dwp_core ), ]


printlog( "Computing time elapsed" )
alldata <- alldata[ order( alldata[ 'EPOCH_TIME' ] ), ]
alldata <- ddply( alldata, .( treatment, rep ), mutate, ELAPSED_TIME=( EPOCH_TIME-EPOCH_TIME[ 1 ] )/60, .progress="text" )


fd <- read_csv( "core_data.csv", INPUT_DIR )
printlog( "Merging Picarro and field data..." )
alldata <- merge( alldata, fd )

print( summary( alldata ) )

alldata$dwp_core <- as.factor( alldata$dwp_core )

p_ch4 <- qplot( ELAPSED_TIME, CH4_dry, data=alldata, color=rep )
p_ch4 <- p_ch4 + facet_grid( treatment~., scales="free" ) + scale_color_discrete( "DWP core" )
print( p_ch4 )
saveplot( "summary_ch4_allreps" )

p_co2 <- qplot( ELAPSED_TIME, CO2_dry, data=alldata, color=rep )
p_co2 <- p_co2 + facet_grid( treatment~., scales="free" ) + scale_color_discrete( "DWP core" )
print( p_co2 )
saveplot( "summary_co2_allreps" )


	printlog( "Focusing on rep 1 only..." )
	alldata1 <- alldata[ alldata$rep=="Rep 1", ]


p_ch4r1 <- qplot( ELAPSED_TIME, CH4_dry, data=alldata1, geom="line", group=1, size=I( 2 ), color=dwp_core )
p_ch4r1 <- p_ch4r1 + facet_grid( treatment~., scales="free" )
print( p_ch4r1 )
saveplot( "summary_ch4_rep1" )

p_co2r1 <- qplot( ELAPSED_TIME, CO2_dry, data=alldata1, geom="line", group=1, size=I( 2 ), color=dwp_core )
p_co2r1 <- p_co2r1 + facet_grid( treatment~., scales="free" ) 
print( p_co2r1 )
saveplot( "summary_co2_rep1" )


saveplot( "summary_ch4_rep1_closeup", p_ch4r1 + xlim( c( 1000,1080 ) ) )

#print( p_co2r1 + xlim( c( 1000,1080 ) ) )
#saveplot( "summary_co2_rep1_closeup" )

stop('ok')


printlog( "Computing fluxes..." )
fluxes <- ddply( alldata, .( filename, Plot ), .fun=compute_flux )

print( summary( fluxes ) )

p <- ggplot( fluxes, aes( Day, flux, group=Plot, colour=Plot ) ) + geom_point() + geom_line()
print( p )
saveplot( "flux_summary" )

printlog( SEPARATOR )
printlog( "Saving flux data..." )
savedata( alldata )
printlog( "Saving flux data..." )
savedata( fluxes )

printlog( "All done with", SCRIPTNAME )
print( sessionInfo() )
sink()
