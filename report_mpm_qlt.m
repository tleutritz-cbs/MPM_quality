function report_mpm_qlt(ResPath)
cpn = fullfile(ResPath,'Supplementary');
contr = {'MTw' 'PDw' 'T1w'};
ref.MTw = 5.209;
ref.PDw = 4.796;
ref.T1w = 4.650;
coreg = {'MT2PD' 'T12PD'};
lim.trans = 0.3778;
lim.rot = 0.1137;
% prepare report file
fid = fopen(fullfile(cpn,'quality_report.htm'),'w');
str = ['<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"' ...
    '"http://www.w3.org/TR/html4/strict.dtd">\n' ...
    '<html>\n' ...
    '<head>\n' ...
    ' <title>MPM quality check results for folder %s</title>\n' ...
    '</head>\n' ...
    '<h1>Quality check results for folder %s</h1>\n'
    ];
fprintf(fid,str,ResPath,ResPath);
fprintf(fid,['<h2>Table of contents</h2>\n' ...
    '<a href="#qa_meas">1. hMRI-toolbox QA</a><br>\n' ...
    '<a href="#mpm_qlt">2. MPM quality</a><br>\n' ...
    '<a href="#med_qlt">3. Multi-echo data quality</a>\n' ...
    '<h2><a id="qa_meas">1. hMRI-toolbox QA</a></h2>\n']);
tmp = get_metadata(fullfile(cpn,'hMRI_map_creation_quality_assessment.json'));
tmp = tmp{1,1};
WMPDerr = 100*tmp.PD.SD/tmp.PD.mean;
fprintf(fid,'Error estimate within PD WM (measure for the accuracy of the bias-field correction): %5.2f%%<br>\n',WMPDerr);
fprintf(fid,'For comparison: within THS at ZH the average value was 10.25%%<br>\n');
if WMPDerr > 10.25
    fprintf(fid,'<font color="red">The actual value is higher by %5.2f p.u.</font><br>\n',WMPDerr-10.25);
end
fprintf(fid,['<p>Standard deviations within white matter of the R2* maps ' ...
    'calculated from each individual contrast (measure for intra-scan motion):\n<ul>\n']);
for cc = contr
    cn =char(cc);
    cv = tmp.SDR2s.(cn);
    dev = 100*(cv-ref.(cn))/ref.(cn);
    if cv > ref.(cn)
        fprintf(fid,'<li><font color="red">from %s: %5.2f (deviation of %5.2f%% to THS)</font></li>\n',cn,cv,dev);
    else
        fprintf(fid,'<li>from %s: %5.2f (deviation of %5.2f%% to THS)</li>\n',cn,cv,dev);
    end
end
fprintf(fid,'</ul></p><p>Inter-scan motion (parameters for coregistration of contrasts to PDw):\n<ul>\n');
for cg = coreg
    cn = char(cg);
    cor = tmp.ContrastCoreg.(cn);
    trans = cor(1:3);
    rot = cor(4:6);
    dev = 100*(rms(trans)-lim.trans)/lim.trans;
    if rms(trans) > lim.trans
        fprintf(fid,'<li><font color="red">%s translation (rms): %5.2f mm (deviation of %5.2f%% to THS)</font></li>\n',cn,rms(trans),dev);
    else
        fprintf(fid,'<li>%s translation (rms): %5.2f mm (deviation of %5.2f%% to THS)</li>\n',cn,rms(trans),dev);
    end
    dev = 100*(rad2deg(rms(rot))-lim.rot)/lim.rot;
    if rad2deg(rms(rot)) > lim.rot
        fprintf(fid,'<li><font color="red">%s rotation (rms): %5.2f ° (deviation of %5.2f%% to THS)</font></li>\n',cn,rad2deg(rms(rot)),dev);
    else
        fprintf(fid,'<li>%s rotation (rms): %5.2f ° (deviation of %5.2f%% to THS)</li>\n',cn,rad2deg(rms(rot)),dev);
    end
end
fprintf(fid,'</ul></p>\n<h2><a id="mpm_qlt">2. MPM quality</a></h2>\n');
fprintf(fid,'<p><img width=800 src="%s" alt="MPM quality 4ROIs"><br>\n',fullfile(ResPath,'MPM_Quality4','MPMqlt_full.png'));
fprintf(fid,['Figure 1: Mean values (top row) and Coefficient of Variation ' ...
    '(CoV; bottom row) for longitudinal and effective transverse relaxation R1 and R2* (in 1/s), ' ...
    'Magnetization Transfer saturation MTsat (in p.u.), proton density PD (in p.u.), and T1-weighted images '...
    'within brain grey matter (GM), caudate nucleus (CN), white matter (WM), and corpus callosum (CC).</p>\n']);
fprintf(fid,'<h2><a id="med_qlt">3. Multi-echo data quality</a></h2>\n(Will be added later.)\n');
fprintf(fid,'</body></html>');
end
