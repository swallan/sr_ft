load.if.needed <-
function(filenames, do.print = T)
{
	# You may wish to insert
	# load.if.needed("qprobs.o")
	# into the definitions of qtukey & ptukey
	# Otherwise use dyn.load("qprobs.o") as required
	if(!exists(".dyn.loaded", frame = 0))
		assign(".dyn.loaded", 	#, mode = "character"))
		"dyn loaded files:", frame = 0)
	if(!exists(".dyn.symbols"))
		assign(".dyn.symbols", 	# , mode = "character"))
		"dyn loaded symbols:", frame = 0)
	ind <- match(.dyn.loaded, filenames, nomatch = 0)
	ind <- ind[ind != 0]
	if(length(ind) == 0)
		loaded <- filenames
	else loaded <- filenames[ - ind]
	if(length(loaded) > 0 && do.print == T)
		print(loaded)
	symbols <- dyn.load(loaded)
	assign(".dyn.loaded", c(.dyn.loaded, loaded), frame = 0)
	assign(".dyn.symbols", c(.dyn.symbols, symbols), frame = 0)
}
