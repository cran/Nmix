print.nmix<-function (x, ...) 
{
cat('Nmix analysis of',x$tag,'\n')
cat('Posterior on k:\n')
np<-length(x$post)+min(0,3-match(FALSE,0==rev(x$post)))
cat(x$post[1:np],fill=60)
cat('nsweep:',x$nsweep,'nburnin:',x$nburnin,'seed:',x$seed,
        'moves:',x$moves,'\n')
invisible(x)
}
