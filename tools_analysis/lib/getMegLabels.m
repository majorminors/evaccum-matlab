function labels = getMegLabels(region,hemisphere,sensorType)

% labels by lobe(ish) from midline to lateral edge, anterior to posterior
% so e.g. .megmags.frontal.left(1:4) will get you most anterior midline
% magneto sensor followed up second most anterior midline sensor, followed
% by third, then first of the next lateral row of sensors

validValues = {'temporal', 'parietal', 'frontal', 'occipital'};
if nargin < 1 || isempty(region)
    error('need a region')
    return
else
    if ~ismember(region, validValues)
        errorMsg = ['region must be one of:'];
        for i = 1:numel(validValues)
            errorMsg = [errorMsg ' ' validValues{i}];
        end
        error(errorMsg);
        return
    end
end
clear validValues errorMsg

validValues = {'left' 'right'};
if nargin < 2 || isempty(hemisphere)
    hemisphere = validValues;
else
    if ~ismember(region, validValues)
        errorMsg = ['hemisphere must be one of:'];
        for i = 1:numel(validValues)
            errorMsg = [errorMsg ' ' validValues{i}];
        end
        error(errorMsg);
        return
    end
end
clear validValues errorMsg

validValues = {'megmags' 'megplanars'};
if nargin < 3 || isempty(sensorType)
    sensorType = validValues;
else
    if ~ismember(sensorType, validValues)
        errorMsg = ['sensorType must be one of:'];
        for i = 1:numel(validValues)
            errorMsg = [errorMsg ' ' validValues{i}];
        end
        error(errorMsg);
        return
    end
end
clear validValues errorMsg

sensors{1,1} = {'occipital'};
sensors{1,2} = {'right'};
sensors{1,3} = {'megplanars'};
sensors{1,4} = {'MEG2032' 'MEG2033' 'MEG2122' 'MEG2123' 'MEG2132' 'MEG2133' 'MEG2342' 'MEG2343' 'MEG2332' 'MEG2333' 'MEG2312' 'MEG2313' 'MEG2322' 'MEG2323' 'MEG2512' 'MEG2513' 'MEG2542' 'MEG2543' 'MEG2432' 'MEG2433' 'MEG2522' 'MEG2523' 'MEG2532' 'MEG2533'};
sensors{2,1} = {'occipital'};
sensors{2,2} = {'right'};
sensors{2,3} = {'megmags'};
sensors{2,4} = {'MEG2031' 'MEG2121' 'MEG2131'  'MEG2341' 'MEG 2331' 'MEG2311' 'MEG2321' 'MEG2511' 'MEG2541' 'MEG2431' 'MEG2521' 'MEG2531'};

sensors{3,1} = {'occipital'};
sensors{3,2} = {'left'};
sensors{3,3} = {'megplanars'};
sensors{3,4} = {'MEG2042' 'MEG2043' 'MEG2112' 'MEG2113' 'MEG2142' 'MEG2143' 'MEG1912' 'MEG1913' 'MEG1922' 'MEG1923' 'MEG1932' 'MEG1933' 'MEG1942' 'MEG1943' 'MEG1732' 'MEG1733' 'MEG1742' 'MEG1743' 'MEG1642' 'MEG1643' 'MEG1722' 'MEG1723' 'MEG1712' 'MEG1713'};
sensors{4,1} = {'occipital'};
sensors{4,2} = {'left'};
sensors{4,3} = {'megmags'};
sensors{4,4} = {'MEG2041' 'MEG2111' 'MEG2141' 'MEG1911' 'MEG1921' 'MEG1931' 'MEG1941' 'MEG1731' 'MEG1741' 'MEG1641' 'MEG1721' 'MEG1711'};

sensors{5,1} = {'parietal'};
sensors{5,2} = {'right'};
sensors{5,3} = {'megplanars'};
sensors{5,4} = {'MEG1042' 'MEG1043' 'MEG0722' 'MEG0723' 'MEG0732' 'MEG0733' 'MEG2242' 'MEG2243' 'MEG2022' 'MEG2033' 'MEG1112' 'MEG1113' 'MEG1142' 'MEG1143' 'MEG2212' 'MEG2213' 'MEG2232' 'MEG2233' 'MEG1122' 'MEG1123' 'MEG1132' 'MEG1133' 'MEG2222' 'MEG2223' 'MEG2442' 'MEG2443'};
sensors{6,1} = {'parietal'};
sensors{6,2} = {'right'};
sensors{6,3} = {'megmags'};
sensors{6,4} = {'MEG1041' 'MEG0721' 'MEG0731' 'MEG2241' 'MEG2021' 'MEG1111' 'MEG1141' 'MEG2211' 'MEG2231' 'MEG1121' 'MEG1131' 'MEG2221'  'MEG2441'};

sensors{7,1} = {'parietal'};
sensors{7,2} = {'left'};
sensors{7,3} = {'megplanars'};
sensors{7,4} = {'MEG0632' 'MEG0633' 'MEG0712' 'MEG0713' 'MEG0742' 'MEG0743' 'MEG1832' 'MEG1833' 'MEG2012' 'MEG2013' 'MEG0422' 'MEG0423' 'MEG0432' 'MEG0433' 'MEG1822' 'MEG1823' 'MEG1842' 'MEG1843' 'MEG0412' 'MEG0413' 'MEG0442' 'MEG0443' 'MEG1812' 'MEG1813' 'MEG1632' 'MEG1633'};
sensors{8,1} = {'parietal'};
sensors{8,2} = {'left'};
sensors{8,3} = {'megmags'};
sensors{8,4} = {'MEG0631' 'MEG0711' 'MEG0741' 'MEG1831' 'MEG2011' 'MEG0421' 'MEG0431' 'MEG1821' 'MEG1841' 'MEG0411' 'MEG0441' 'MEG1811' 'MEG1631'};

sensors{9,1} = {'temporal'};
sensors{9,2} = {'right'};
sensors{9,3} = {'megplanars'};
sensors{9,4} = {'MEG1312' 'MEG1313' 'MEG1342' 'MEG1343' 'MEG2412' 'MEG2413' 'MEG1322' 'MEG1323' 'MEG1332' 'MEG1333' 'MEG2422' 'MEG2423' 'MEG1442' 'MEG1443' 'MEG2612' 'MEG2613' 'MEG2642' 'MEG2643' 'MEG1422' 'MEG1423' 'MEG1432' 'MEG1433' 'MEG2622' 'MEG2623' 'MEG2632' 'MEG2633'};
sensors{10,1} = {'temporal'};
sensors{10,2} = {'right'};
sensors{10,3} = {'megmags'};
sensors{10,4} = {'MEG1311' 'MEG1341' 'MEG2411' 'MEG1321' 'MEG1331' 'MEG2421' 'MEG1441' 'MEG2611' 'MEG2641' 'MEG1421' 'MEG1431' 'MEG2621' 'MEG2631'};

sensors{11,1} = {'temporal'};
sensors{11,2} = {'left'};
sensors{11,3} = {'megplanars'};
sensors{11,4} = {'MEG0222' 'MEG0223' 'MEG0232' 'MEG0233' 'MEG1622' 'MEG1623' 'MEG0212' 'MEG0213' 'MEG0242' 'MEG0243' 'MEG1612' 'MEG1613' 'MEG0132' 'MEG0133' 'MEG1512' 'MEG1513' 'MEG1522' 'MEG1523' 'MEG0112' 'MEG0113' 'MEG0142' 'MEG0143' 'MEG1542' 'MEG1543' 'MEG1532' 'MEG1533'};
sensors{12,1} = {'temporal'};
sensors{12,2} = {'left'};
sensors{12,3} = {'megmags'};
sensors{12,4} = {'MEG0221' 'MEG0231' 'MEG1621' 'MEG0211' 'MEG0241' 'MEG1611' 'MEG0131' 'MEG1511' 'MEG1521' 'MEG0111' 'MEG0141' 'MEG1541' 'MEG1531'};

sensors{13,1} = {'frontal'};
sensors{13,2} = {'right'};
sensors{13,3} = {'megplanars'};
sensors{13,4} = {'MEG0812' 'MEG0813' 'MEG1012' 'MEG1013' 'MEG0912' 'MEG0913' 'MEG0942' 'MEG0943' 'MEG1022' 'MEG1023' 'MEG1032' 'MEG1033' 'MEG0922' 'MEG0923' 'MEG0932' 'MEG0933' 'MEG1242' 'MEG1243' 'MEG1212' 'MEG1213' 'MEG1232' 'MEG1233' 'MEG1222' 'MEG1223' 'MEG1412' 'MEG1413'};
sensors{14,1} = {'frontal'};
sensors{14,2} = {'right'};
sensors{14,3} = {'megmags'};
sensors{14,4} = {'MEG0811' 'MEG1011' 'MEG0911' 'MEG0941' 'MEG1021' 'MEG1031' 'MEG0921' 'MEG0931' 'MEG1241' 'MEG1211' 'MEG1231' 'MEG1221' 'MEG1411'};

sensors{15,1} = {'frontal'};
sensors{15,2} = {'left'};
sensors{15,3} = {'megplanars'};
sensors{15,4} = {'MEG0822' 'MEG0823' 'MEG0622' 'MEG0623' 'MEG0522' 'MEG0523' 'MEG0532' 'MEG0533' 'MEG0612' 'MEG0613' 'MEG0642' 'MEG0643' 'MEG0512' 'MEG0513' 'MEG0542' 'MEG0543' 'MEG0332' 'MEG0333' 'MEG0312' 'MEG0313' 'MEG0322' 'MEG0323' 'MEG0342' 'MEG0343' 'MEG0122' 'MEG0123'};
sensors{16,1} = {'frontal'};
sensors{16,2} = {'left'};
sensors{16,3} = {'megmags'};
sensors{16,4} = {'MEG0821' 'MEG0621' 'MEG0521' 'MEG0531' 'MEG0611' 'MEG0641' 'MEG0511' 'MEG0541' 'MEG0331' 'MEG0311' 'MEG0321' 'MEG0341' 'MEG0121'};

labels = sensors(cellfun(@(x) strcmp(x,region), sensors(:,1)) & cellfun(@(x) ismember(x,hemisphere), sensors(:,2)) & cellfun(@(x) ismember(x,sensorType), sensors(:,3)), 4);
labels = cat(2,labels{:});

end
