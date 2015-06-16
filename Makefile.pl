use Cwd;

$repo_path = 'd:/REPO';
$gitTags_file = cwd.'/gitTags';
$reps = [   {path => $repo_path.'/my_lib', tag_mask => '102001_MY_LIB[\w\d]+'}  ];

open(LABELS,"$gitTags_file");
@strs = <LABELS>;
foreach (@strs)
{
	chomp($_);
    foreach $rep (@$reps)
    {
        if ($_ =~ m/\A($$rep{tag_mask})\Z/)
        {
            if (exists($$rep{tag}) && !($$rep{tag} eq $1))
            {
                print "ERROR: diplicate tag found for $$rep{path} $$rep{tag} -> $1";
                close(LABELS); exit;
            }
            else
            {
                if (!($$rep{tag} eq $1))
                {
                    $$rep{tag} = $1;
                    print "Found label for $$rep{path}: $$rep{tag}\n";                            
                }
            }
        }
    }
}
close(LABELS);

foreach $rep (@$reps)
{
    if (!exists($$rep{tag}))
    {
        print "ERROR: no tag found for $$rep{path}";
        exit;
    }
}
print "Applying tags...\n";
foreach $rep (@$reps)
{
    chdir("$$rep{path}");
    $cmd = "git checkout $$rep{tag}";    
    print "$$rep{path} -> \n";
    system($cmd); 
}
