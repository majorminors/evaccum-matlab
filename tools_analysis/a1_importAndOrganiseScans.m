function a1_importAndOrganiseScans(thisSubject)

baseF{sb} = [fld_tar,num2str(subjs(sb)),'/'];

ID = getMEGID(sprintf('DM_evaccumpilot_%s',num2str(subjs(sb))));
tidyup_evaccum(baseF{sb},ID,overwrite);

return
end