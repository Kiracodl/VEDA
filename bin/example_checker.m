%the problem: I sometimes change the structure in the gui but forget to change the example files.
%I could try to run every example at every release and make sure everything is okay.
%easier to run a  script that scans an example xml file and tries to see if there is any tag that 
%is not also present in the corresponding tool xml file.  

function tags_lister(examplefile, toolfile) 

xdoc = xmlread(examplefile);

tool = xmlread(toolfile);

k = 1;
for nam = { 'number', 'boolean', 'string', 'choice'}
    num = xdoc.getElementsByTagName(nam).getLength;
    for i = 0:num-1
        c = xdoc.getElementsByTagName(nam).item(i);
        id = c.getAttribute('id');
        path1 = walk_up(c);        
        
        list_ex{k} = path1;
        k=k+1;
    end
end

k = 1;
for nam = { 'number', 'boolean', 'string', 'choice'}
    num = tool.getElementsByTagName(nam).getLength;
    for i = 0:num-1
        c = tool.getElementsByTagName(nam).item(i);
        id = c.getAttribute('id');
        path1 = walk_up(c);        
        
        list_tool{k} = path1;
        k=k+1;
    end
end

list_d = setdiff( list_ex, list_tool);
disp(list_d');
end

function s = walk_up( c )
tn = c.getTagName;
id = c.getAttribute('id');
if ( strcmp(tn, 'input') == 0)
    s1 = walk_up(c.getParentNode);
else
    s1 = [];
end
s = [s1 '.' char(tn) '(' char(id) ')' ];
end