function [s,desc]=readhdf5complex(fn)
% [s,desc]=readhdf5complex(fn)
% Reads and HDF5 file, consisting of attributes and datasets with as much
% complexity as is present -- i.e. folders within folders.  Differs from
% readhdf5both which assumes one layer of structure.
%
% Optional output desc holds a list of all the fields as (1) hdf5 names and as
% (2) field names, (3) data size, (4) isstr and (5) isattribute
%
% NOTE: field names with spaces become double underscores, e.g. 
% 'a b' becomes 'a__b' so it can be legal field name in a structure.
% writehdf5complex converts double underscores back to spaces.
% 
% Updated 04/05/2011 to fix the problem of naming fields of the structure
% starting with a number, i.e. 35_ngs_wfs_quad, which when attempting to
% place in a structure will give you an error, because variables or
% fieldnames cannot start with an integer -  must be an alphanumeric
% character.  I put an X in front of the name, so now it will return the
% 35_ngs_wfs_quad info in a field labeled X35_ngs_wfs_quad.  - LAK

global desc ptr currfld

if ~isstruct(fn)
    if isstr(fn)
        fn = hdf5info(fn);
    else
        error('Input must be a structure or a filename')
    end
end
g=fn.GroupHierarchy;
desc=cell(100,5);  % leaving room
ptr=1;
currfld='.';
s=proc(g);
desc=desc(1:ptr-1,:);
return

%%%%%%%%%%%%%%%%%%%%%%%%
function s=proc(g)
%%%%%%%%%%%%%%%%%%%%%%%%
global desc ptr currfld
% ATTRIBUTES
attr=g.Attributes;
if ~isempty(attr)
    for i=1:length(attr)
        a=attr(i);
        name=nameonly(a.Name);
        aa=a.Value;
        if isa(aa,'hdf5.h5string')|isa(aa,'hdf5.h5array')
            aa=aa.Data;
        end
        s.(name) = aa;
        desc(ptr,:)={a.Name [currfld name] size(aa) isstr(aa) 1};
        ptr=ptr+1;
    end
end
% DATASETS
data=g.Datasets;
if ~isempty(data)
    for i=1:length(data)
        d=data(i);
        name = nameonly(d.Name);
        dd=hdf5read(d);
        if isa(dd,'hdf5.h5string')
            dd=dd.Data;
            f=find(dd==char(10));
            if ~isempty(f)
                dd=reshape(dd,f(1),length(dd)/f(1))';
                dd=dd(:,1:end-1); % clip off the char(10)'s
            end
        elseif isa(dd,'hdf5.h5array')
            dd=dd.Data;
        end
        s.(name) = dd;
        desc(ptr,:)={d.Name [currfld name] size(dd) isstr(dd) 0};
        ptr=ptr+1;
    end
end
% GROUPS
grp=g.Groups;
if ~isempty(grp)
    currfld0=currfld;
    for i=1:length(grp)
        g=grp(i);
        name = nameonly(g.Name);
        currfld=[currfld0 name '.'];
        try
            s.(name) = proc( g );
        catch
            name = ['X' name];
            s.(name) = proc( g );
        end
    end
end
try
    isstruct(s);
catch
    s=struct([]);
end
return
        

%%%%%%%%
function n = nameonly(nn)
%%%%%%%%
% FIND "LAST" NAME
f=find(nn=='/');
n=nn((f(end)+1):end);
% REPLACE SPACES WITH DOUBLE UNDERSCORES __
f=find(' '==n);
if ~isempty(f)
    for i=1:length(f)
        f=find(' '==n);
        f0=f(1);
        n=[n(1:f0-1) '__' n(f0+1:end)];
    end
end
% REPLACE NON-ALPHANUMERICS WITH SINGLE UNDERSCORES
f=~isalpha(n);
n(f)='_';
return

%%%%%%%%%%%%
function f = isalpha(str)
%%%%%%%%%%%%
% returns ones for alphanumerics
f = (64<str & str<91)|(96<str & str<123)|(47<str & str<58);
return
