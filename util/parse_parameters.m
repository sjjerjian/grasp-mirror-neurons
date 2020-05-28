function p=parse_parameters(v,var_names,defaults)
%parses parameters given as a sequence of pairs of variable name and its value
%assigns defaults to the variables without given values

for i=1:length(var_names),
   p.(var_names{i})=defaults{i};
end
if floor(length(v)/2)~=length(v)/2, v(end)=[];end
p.nonrecognized={};
for i=1:2:length(v),
   s=v{i};
   if ischar(s), 
       t=find(strcmpi(s,var_names));
%       t=strmatch(lower(s),var_names,'exact')
      if ~isempty(t),
         p.(var_names{t})=v{i+1};
      else
          p.nonrecognized{end+1}=s;
      end
   else
      fprintf('parse_parameters:Something is wrong with optional arguments, just ignoring them\n');
   end
end
end