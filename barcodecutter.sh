awk '
BEGIN{
    x="'$1'"
    match(x, /^(N+)([ACGTU]+)(N+)$/, bc);
    l=length(x);
    r=bc[1, "length"];
    tag=bc[2]
    print "l="l" r="r" tag="tag
}
{
    if( match($0, "^.{"r"}"tag) ) {
        print last;
        print substr($0, l);
        getline;
        print;
        getline;
        print substr($0, l)
    } else{
        last=$0
    }
}'
