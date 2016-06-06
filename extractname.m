function [expression,FieldNames] = extractname(extension)

fprintf(['This module enables the segmentation of the filenames of the\n' ...
    'seperate images. An example could be:\n'])

fprintf(['\nFilename: generalinfo_T014 \n'...
    'Where ''T014'' is timepoint 14 and we want to exclude the ''generalinfo_'' part\n' ...
    '\nThis module would work as follows:\n\n' ...
    'Is there a defining character? (y/n)\n' ...
    'y\n\n' ...
    'What is the defining character?\n' ...
    'T\n\n' ...
    'Is the indicator a number? (y/n)\n' ...
    'y\n\n' ...
    'With what parameter do you want to save this indicator?(e.g. t or RowLabel)\n' ...
    'time\n\n' ...
    'The output of this example would be:\n\n' ...
    'expression = ^.*T(\\d+).tiff\n' ...
    'storage = ''time''}\n\n' ...
    'Now press any key to start the module\n\n'])
    pause

FieldNames = {};

%% SPECIFY FIRST PART STRING

promt = 'Does the filename start with an unrelevant part? (y/n)\n';
answer = input(promt,'s')';

if strcmp(answer,'y')
    promt = ['\nWhat are the last characters of the unrelevant part?\n' ...
        '- Can not be present in earlier part of the string\n' ...
        '- Is sensitive to lower/upper case\n'];
    answer = input(promt,'s')';
    
    answer = fixanswer(answer);
    
    expression = ['^.*' answer];
else
    expression = '';
end

%% LOOP THROUGH THE NEXT PARTS OF THE STRING

while (1)
    
    % Is there a character that is not included but indicates information?
    promt = '\nIs there a defining character? (y/n)\n';
    answer = input(promt,'s')';
    if strcmp(answer,'y')
        promt = '\nWhat is the defining character?\n';
        answer = input(promt,'s')';
        answer = fixanswer(answer);
        expression = [expression answer];
    end
    
    % What type of information are we looking for?
    promt = '\nIs the indicator a number? (y/n)\n';
    answer = input(promt,'s')';
    if strcmp(answer,'y')
        expression = [expression '(\d+)'];
    else
        promt = '\nWhat char are we accepting? (e.g. a-b or A-H)\n';
        answer = input(promt,'s')';
        answer = fixanswer(answer);
        expression = [expression '([' answer '])'];
    end
    
    % What information does this indicate?
    promt = ['\nWith what parameter do you want to save this indicator?' ...
        '(e.g. t or RowLabel)\n'];
    answer = input(promt,'s')';
    answer = fixanswer(answer);
    FieldNames{end+1,1} = answer;
    
    % Do you want to continue?
    promt = '\nDo you want to add another character? (y/n)\n';
    answer = input(promt,'s')';
    if strcmp(answer,'n')
        break
    end
    
end

%% ADD FILE EXTENSION

expression = [expression '.' extension];


end