classdef GeoRank_3D_App < matlab.apps.AppBase
%GEORANK_3D_APP  GeoRank_3D interactive GUI — complete validation workflow.
%
%   LAUNCH:   GeoRank_3D_App
%   INSTALL:  MATLAB Apps tab > Install App > select this file
%
%   Implements all 7 validation sources from the original documentation:
%     1. Point Cloud  (GeoRank_3D pipeline)
%     2. Manual Sets  (field compass CSV)
%     3. qFacet FM    (Fast Marching export CSV)
%     4. qFacet KD    (KD-tree export CSV)
%     5. DSE Sets     (Discontinuity Set Extractor CSV)
%     6. Observed Blocks (field block measurements CSV)
%     7. RocFall3 Output (Rocscience RocFall3 CSV)

    %% ---- UI Components -----------------------------------------------
    properties (Access = public)
        UIFigure        matlab.ui.Figure
        LeftPanel       matlab.ui.container.Panel
        RightPanel      matlab.ui.container.Panel

        % Files
        FilesPanel      matlab.ui.container.Panel
        InputField      matlab.ui.control.EditField
        BrowseInBtn     matlab.ui.control.Button
        OutputField     matlab.ui.control.EditField
        BrowseOutBtn    matlab.ui.control.Button

        % Slope params
        SlopePanel      matlab.ui.container.Panel
        SlopeDipSpin    matlab.ui.control.Spinner
        SlopeDDirSpin   matlab.ui.control.Spinner
        FrictionSpin    matlab.ui.control.Spinner

        % Algorithm params
        AlgoPanel       matlab.ui.container.Panel
        MaxPtsFld       matlab.ui.control.NumericEditField
        VoxelFld        matlab.ui.control.NumericEditField
        NormalKFld      matlab.ui.control.NumericEditField
        AngleEpsFld     matlab.ui.control.NumericEditField
        MinClFld        matlab.ui.control.NumericEditField
        PlanarMinFld    matlab.ui.control.NumericEditField
        MCRunsFld       matlab.ui.control.NumericEditField
        PlanarWinFld    matlab.ui.control.NumericEditField
        ToppWinFld      matlab.ui.control.NumericEditField

        % Validation sources panel
        ValSrcPanel     matlab.ui.container.Panel
        ManualFld       matlab.ui.control.EditField
        ManualBtn       matlab.ui.control.Button
        QfFMFld         matlab.ui.control.EditField
        QfFMBtn         matlab.ui.control.Button
        QfKDFld         matlab.ui.control.EditField
        QfKDBtn         matlab.ui.control.Button
        DSEFld          matlab.ui.control.EditField
        DSEBtn          matlab.ui.control.Button
        ObsFld          matlab.ui.control.EditField
        ObsBtn          matlab.ui.control.Button
        RF3Fld          matlab.ui.control.EditField
        RF3Btn          matlab.ui.control.Button
        RunValBtn       matlab.ui.control.Button

        % Run / status
        RunBtn          matlab.ui.control.Button
        StatusLbl       matlab.ui.control.Label

        % Results tabs
        TabGroup        matlab.ui.container.TabGroup
        Tab_Sets        matlab.ui.container.Tab
        Tbl_Sets        matlab.ui.control.Table
        Tab_Pole        matlab.ui.container.Tab
        Ax_Pole         matlab.ui.control.UIAxes
        Tab_Rel         matlab.ui.container.Tab
        Ax_Rel          matlab.ui.control.UIAxes
        Tab_Wedge       matlab.ui.container.Tab
        Tbl_Wedge       matlab.ui.control.Table
        Tab_Disp        matlab.ui.container.Tab
        Tbl_DispSets    matlab.ui.control.Table
        Tbl_DispWedge   matlab.ui.control.Table
        Tab_MultiVal    matlab.ui.container.Tab
        Tbl_MultiVal    matlab.ui.control.Table
        Ax_DevHeat      matlab.ui.control.UIAxes
        Tab_BlockComp   matlab.ui.container.Tab
        Tbl_BlockComp   matlab.ui.control.Table
        Ax_BlockHist    matlab.ui.control.UIAxes
        Tab_RocFall     matlab.ui.container.Tab
        Tbl_RocStats    matlab.ui.control.Table
        Ax_RocHist      matlab.ui.control.UIAxes
        Tab_Log         matlab.ui.container.Tab
        LogArea         matlab.ui.control.TextArea
    end

    properties (Access = private)
        Res_
    end

    %% ---- Callbacks ---------------------------------------------------
    methods (Access = private)

        function browseFile(app, field, filter, title)
            [f,p] = uigetfile(filter, title);
            if ~isequal(f,0); field.Value = fullfile(p,f); end
        end

        function onBrowseIn(app,~)
            app.browseFile(app.InputField,...
                {'*.xyz;*.txt;*.csv;*.ply','Point Cloud'},'Select point cloud');
        end
        function onBrowseOut(app,~)
            p=uigetdir(pwd,'Select output folder');
            if ~isequal(p,0); app.OutputField.Value=p; end
        end
        function onBrowseManual(app,~)
            app.browseFile(app.ManualFld,{'*.csv;*.txt','CSV'},'Manual field compass sets CSV');
        end
        function onBrowseQfFM(app,~)
            app.browseFile(app.QfFMFld,{'*.csv;*.txt','CSV'},'qFacet Fast Marching export CSV');
        end
        function onBrowseQfKD(app,~)
            app.browseFile(app.QfKDFld,{'*.csv;*.txt','CSV'},'qFacet KD-tree export CSV');
        end
        function onBrowseDSE(app,~)
            app.browseFile(app.DSEFld,{'*.csv;*.txt','CSV'},'DSE export CSV');
        end
        function onBrowseObs(app,~)
            app.browseFile(app.ObsFld,{'*.csv;*.txt','CSV'},'Observed blocks CSV');
        end
        function onBrowseRF3(app,~)
            app.browseFile(app.RF3Fld,{'*.csv;*.txt','CSV'},'RocFall3 output CSV');
        end

        function onRun(app,~)
            inFile = strtrim(app.InputField.Value);
            outDir = strtrim(app.OutputField.Value);
            if isempty(inFile)||~isfile(inFile)
                uialert(app.UIFigure,'Select a valid point cloud file.','No Input'); return;
            end
            if isempty(outDir); outDir=fullfile(fileparts(inFile),'GeoRank_3D_output'); end

            opts = GeoRank_3D_default_options();
            opts.maxPoints      = max(1, round(app.MaxPtsFld.Value*1e6));
            opts.slopeDip       = app.SlopeDipSpin.Value;
            opts.slopeDipDir    = app.SlopeDDirSpin.Value;
            opts.frictionAngle  = app.FrictionSpin.Value;
            opts.voxelSize      = app.VoxelFld.Value;
            opts.normalK        = round(app.NormalKFld.Value);
            opts.angularEpsDeg  = app.AngleEpsFld.Value;
            opts.minClusterPts  = round(app.MinClFld.Value);
            opts.planarityMin   = app.PlanarMinFld.Value;
            opts.monteCarloRuns = round(app.MCRunsFld.Value);
            opts.planarKinWin   = app.PlanarWinFld.Value;
            opts.topplingKinWin = app.ToppWinFld.Value;

            app.setStatus('Running GeoRank_3D ...','#B06000');
            app.RunBtn.Enable='off';
            app.log(repmat('-',1,50));
            app.log(sprintf('Run: %s', datestr(now,'yyyy-mm-dd HH:MM:SS')));
            app.log(sprintf('Input:  %s', inFile));
            drawnow;

            try
                res = GeoRank_3D(inFile, outDir, opts);
                app.Res_ = res;
                app.showResults(res);
                app.setStatus(sprintf('Done — %d sets extracted.',height(res.setTable)),'#005500');
                app.log(sprintf('Complete: %d sets | %d patches | %d wedge pairs',...
                    height(res.setTable),height(res.patchTable),height(res.wedgeTable)));
                if isfield(res,'setDispTable') && ~isempty(res.setDispTable)
                    app.log('Dispersion-aware kinematics (population vs mean):');
                    for kk=1:height(res.setDispTable)
                        app.log(sprintf('  Set %d: planar %.1f%% / toppling %.1f%% of population feasible (mean %d/%d)',...
                            res.setDispTable.SetID(kk),100*res.setDispTable.PlanarFeasFrac(kk),...
                            100*res.setDispTable.TopplingFeasFrac(kk),...
                            res.setDispTable.PlanarFeasMean(kk),res.setDispTable.TopplingFeasMean(kk)));
                    end
                    if isfield(res,'wedgeDispTable') && ~isempty(res.wedgeDispTable)
                        for kk=1:height(res.wedgeDispTable)
                            app.log(sprintf('  Wedge %d-%d: %.1f%% of intersections feasible (mean %d)',...
                                res.wedgeDispTable.SetA(kk),res.wedgeDispTable.SetB(kk),...
                                100*res.wedgeDispTable.WedgeFeasFrac(kk),res.wedgeDispTable.WedgeFeasMean(kk)));
                        end
                    end
                end
            catch ME
                app.setStatus(['ERROR: ' ME.message],'#AA0000');
                app.log(['ERROR: ' ME.message]);
                uialert(app.UIFigure,ME.message,'Analysis Error');
            end
            app.RunBtn.Enable='on';
        end

        function onRunValidation(app,~)
            if isempty(app.Res_)||isempty(app.Res_.setTable)
                uialert(app.UIFigure,'Run GeoRank_3D analysis first.','No Results'); return;
            end
            res = app.Res_;
            app.setStatus('Running validation ...','#B06000');
            app.RunValBtn.Enable='off';
            app.log('--- Validation started ---');
            drawnow;

            try
                %-- Build sources struct for set comparison (sources 2-5)
                sources = struct();
                srcDefs = {app.ManualFld.Value,'manual','Manual';
                           app.QfFMFld.Value, 'qfacet_fm','qFacet FM';
                           app.QfKDFld.Value, 'qfacet_kd','qFacet KD';
                           app.DSEFld.Value,  'dse','DSE'};

                anySrc = false;
                for si=1:size(srcDefs,1)
                    fpath=strtrim(srcDefs{si,1});
                    fname=srcDefs{si,2};
                    if ~isempty(fpath)&&isfile(fpath)
                        [dps,dds,lbs,~]=GeoRank_3D_read_reference_csv(fpath,fname);
                        if ~isempty(dps)
                            sources.(fname).dips=dps; sources.(fname).dipDirs=dds; sources.(fname).labels=lbs;
                            app.log(sprintf('Loaded %s: %d sets from %s',srcDefs{si,3},numel(dps),fpath));
                            anySrc=true;
                        end
                    end
                end

                if anySrc
                    [compT, devMat] = GeoRank_3D_validate_multi_source(res.setTable, sources);
                    app.Tbl_MultiVal.Data = compT;
                    app.drawDevHeatmap(devMat, res.setTable, sources);
                    writetable(compT, fullfile(res.outDir,'GeoRank_3D_multi_source_comparison.csv'));
                    app.TabGroup.SelectedTab = app.Tab_MultiVal;
                    app.log('Multi-source comparison written to CSV.');
                else
                    app.log('No set comparison sources loaded (Manual/qFacet FM/qFacet KD/DSE).');
                end

                %-- Block volume comparison (source 6)
                obsFpath = strtrim(app.ObsFld.Value);
                if ~isempty(obsFpath)&&isfile(obsFpath)
                    blockCompT = GeoRank_3D_compare_block_volumes(res.patchTable, obsFpath);
                    if ~isempty(blockCompT)
                        app.Tbl_BlockComp.Data = blockCompT;
                        app.drawBlockComparison(res.patchTable, obsFpath);
                        writetable(blockCompT,fullfile(res.outDir,'GeoRank_3D_block_comparison.csv'));
                        app.log('Block volume comparison written to CSV.');
                    end
                else
                    app.log('Observed blocks file not loaded (source 6).');
                end

                %-- RocFall3 comparison (source 7)
                rf3Path = strtrim(app.RF3Fld.Value);
                if ~isempty(rf3Path)&&isfile(rf3Path)
                    [rocStats, rocT] = GeoRank_3D_parse_rocfall3(rf3Path, res.patchTable);
                    if ~isempty(rocStats)
                        app.showRocStats(rocStats, rocT);
                        writetable(rocT,fullfile(res.outDir,'GeoRank_3D_rocfall3_parsed.csv'));
                        app.log('RocFall3 trajectory stats written to CSV.');
                    end
                else
                    app.log('RocFall3 output file not loaded (source 7).');
                end

                app.setStatus('Validation complete.','#005500');
            catch ME
                app.setStatus(['Validation ERROR: ' ME.message],'#AA0000');
                app.log(['Validation ERROR: ' ME.message]);
            end
            app.RunValBtn.Enable='on';
        end

        %-- Result display helpers ------------------------------------
        function showResults(app,res)
            if ~isempty(res.setTable)
                cols={'SetID','MeanDip_deg','MeanDipDirection_deg','MeanPlanarity',...
                      'PersistenceProxy_m','PlanarFeasible','TopplingFeasible',...
                      'KinematicScore','GeoMechanicalReliability'};
                cols=intersect(cols,res.setTable.Properties.VariableNames,'stable');
                app.Tbl_Sets.Data=res.setTable(:,cols);
            end
            if ~isempty(res.wedgeTable); app.Tbl_Wedge.Data=res.wedgeTable; end
            if isfield(res,'setDispTable') && ~isempty(res.setDispTable); app.Tbl_DispSets.Data=res.setDispTable; end
            if isfield(res,'wedgeDispTable') && ~isempty(res.wedgeDispTable); app.Tbl_DispWedge.Data=res.wedgeDispTable; end
            app.drawPole(res); app.drawRel(res);
            app.TabGroup.SelectedTab=app.Tab_Sets;
        end

        function drawPole(app,res)
            ax=app.Ax_Pole; cla(ax,'reset'); hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');
            if isempty(res.setTable); return; end
            try
                d=res.setTable.MeanDip_deg; a=res.setTable.MeanDipDirection_deg;
                r=sqrt(2).*sind(d./2); xc=r.*sind(a); yc=r.*cosd(a);
                scatter(ax,xc,yc,120,res.setTable.SetID,'filled','MarkerEdgeColor','k','LineWidth',0.6);
                th=linspace(0,2*pi,361);
                plot(ax,sin(th),cos(th),'k-','LineWidth',1.3);
                for rd=[30,60]; rr=sqrt(2)*sind(rd/2);
                    plot(ax,rr*sin(th),rr*cos(th),'--','Color',[.72 .72 .72],'LineWidth',0.5); end
                cds={'N','E','S','W'}; cx=[0,1,0,-1]; cy=[1,0,-1,0];
                for t=1:4
                    plot(ax,[0.95*cx(t),1.06*cx(t)],[0.95*cy(t),1.06*cy(t)],'k-','LineWidth',1.1);
                    text(ax,1.17*cx(t),1.17*cy(t),cds{t},'HorizontalAlignment','center','FontSize',11,'FontWeight','bold');
                end
                for ii=1:height(res.setTable)
                    text(ax,xc(ii)+0.06,yc(ii)+0.06,sprintf('S%d',res.setTable.SetID(ii)),'FontSize',9,'FontWeight','bold');
                end
                colormap(ax,turbo); title(ax,'Schmidt Equal-Area Pole Plot','FontSize',10);
            catch; end; hold(ax,'off');
        end

        function drawRel(app,res)
            ax=app.Ax_Rel; cla(ax,'reset');
            if isempty(res.setTable)||~ismember('GeoMechanicalReliability',res.setTable.Properties.VariableNames); return; end
            try
                hold(ax,'on');
                bar(ax,categorical(string(res.setTable.SetID)),res.setTable.GeoMechanicalReliability,0.6,'FaceColor',[0.18 0.46 0.71]);
                yline(ax,0.5,'--r','LineWidth',1.2);
                text(ax,0.55,0.53,'Threshold=0.5','Units','data','FontSize',8,'Color','r');
                hold(ax,'off');
                xlabel(ax,'Set ID'); ylabel(ax,'Reliability Index'); ylim(ax,[0,1.1]); grid(ax,'on');
                title(ax,'Geomechanical Reliability','FontSize',10);
            catch; end
        end

        function drawDevHeatmap(app, devMat, setTable, sources)
            ax=app.Ax_DevHeat; cla(ax,'reset');
            try
                if isempty(devMat); return; end
                srcNames={'Manual','qFacet FM','qFacet KD','DSE'};
                validSrc=find(~all(isnan(devMat),1));
                if isempty(validSrc); return; end
                D=devMat(:,validSrc); S=srcNames(validSrc);
                imagesc(ax,D); colormap(ax,flipud(summer));
                cb=colorbar(ax); cb.Label.String='Angular deviation [deg]';
                clim(ax,[0,30]);
                set(ax,'XTick',1:numel(S),'XTickLabel',S,'YTick',1:height(setTable),...
                    'YTickLabel',arrayfun(@(x)sprintf('S%d',x),setTable.SetID,'UniformOutput',false));
                title(ax,'Angular deviation heatmap [deg]  (green = pass <10 deg)','FontSize',9);
                for ii=1:size(D,1); for jj=1:size(D,2)
                    if isfinite(D(ii,jj))
                        text(ax,jj,ii,sprintf('%.1f',D(ii,jj)),'HorizontalAlignment','center','FontSize',8,...
                            'Color',ternary(D(ii,jj)<=10,'k','w'));
                    end
                end; end
            catch; end
        end

        function drawBlockComparison(app, patchTable, obsFile)
            ax=app.Ax_BlockHist; cla(ax,'reset'); hold(ax,'on');
            try
                a=patchTable.ProjectedArea_m2; grV=0.20.*sqrt(max(a,0)).^3;
                grV=grV(isfinite(grV)&grV>0);
                obsT=readtable(obsFile,'VariableNamingRule','preserve');
                colN=lower(obsT.Properties.VariableNames);
                vc=find(contains(colN,'vol'),1);
                if isempty(vc); return; end
                obV=double(table2array(obsT(:,vc))); obV=obV(isfinite(obV)&obV>0);
                edges=linspace(min(log10([grV;obV]))-0.5, max(log10([grV;obV]))+0.5, 16);
                histogram(ax,log10(grV),edges,'FaceColor',[0.18 0.46 0.71],'FaceAlpha',0.6,'EdgeColor','none');
                histogram(ax,log10(obV),edges,'FaceColor',[0.85 0.33 0.10],'FaceAlpha',0.6,'EdgeColor','none');
                legend(ax,'GeoRank_3D patches','Observed blocks','Location','northwest');
                xlabel(ax,'log_{10}(Volume [m^3])'); ylabel(ax,'Count');
                title(ax,'Block volume comparison: GeoRank vs Field','FontSize',9);
                grid(ax,'on');
            catch; end; hold(ax,'off');
        end

        function showRocStats(app, rocStats, rocT)
            flds=fieldnames(rocStats);
            rows=cell(numel(flds),2);
            for i=1:numel(flds); rows{i,1}=flds{i}; rows{i,2}=rocStats.(flds{i}); end
            app.Tbl_RocStats.Data=cell2table(rows,'VariableNames',{'Statistic','Value'});
            ax=app.Ax_RocHist; cla(ax,'reset');
            colN=lower(rocT.Properties.VariableNames);
            dc=find(contains(colN,'dist')|contains(colN,'runout'),1);
            if ~isempty(dc)
                d=double(table2array(rocT(:,dc))); d=d(isfinite(d)&d>0);
                if ~isempty(d)
                    histogram(ax,d,14,'FaceColor',[0.85 0.33 0.10],'EdgeColor','w');
                    xlabel(ax,'Runout distance [m]'); ylabel(ax,'Count');
                    title(ax,'RocFall3 runout distribution','FontSize',9); grid(ax,'on');
                end
            end
        end

        function setStatus(app,msg,hexColor)
            app.StatusLbl.Text=msg;
            try
                h=strrep(hexColor,'#','');
                app.StatusLbl.FontColor=double([hex2dec(h(1:2)),hex2dec(h(3:4)),hex2dec(h(5:6))])./255;
            catch; end; drawnow;
        end

        function log(app,msg)
            v=app.LogArea.Value; if ischar(v); v={v}; end
            app.LogArea.Value=[v;{msg}]; drawnow;
        end

    end

    %% ---- Build UI ---------------------------------------------------
    methods (Access = private)
        function build(app)
            W=1400; H=900; LW=360; RW=W-LW-18;

            app.UIFigure=uifigure('Visible','off','Position',[40 40 W H],...
                'Name','GeoRank_3D  —  Reliability-Aware Discontinuity Analysis & Validation',...
                'Color',[0.95 0.95 0.97]);

            %=============================================================
            % LEFT PANEL
            %=============================================================
            app.LeftPanel=uipanel(app.UIFigure,'Title','','BackgroundColor',[0.92 0.93 0.96],...
                'Position',[4 4 LW H-8]);

            uilabel(app.LeftPanel,'Text','GeoRank_3D','FontSize',20,...
                'FontWeight','bold','FontColor',[0.10 0.28 0.60],'Position',[10 H-76 LW-16 40]);
            uilabel(app.LeftPanel,'Text','Reliability-Aware Discontinuity Analysis',...
                'FontSize',9,'FontColor',[0.42 0.42 0.42],'Position',[12 H-98 LW-16 20]);

            y=H-118;

            % Files
            app.FilesPanel=uipanel(app.LeftPanel,'Title','Input / Output','FontWeight','bold','FontSize',10,'Position',[6 y-98 LW-12 102]);
            uilabel(app.FilesPanel,'Text','Input file:','Position',[8 62 80 20],'FontSize',9);
            app.InputField=uieditfield(app.FilesPanel,'text','Position',[90 62 198 22],'FontSize',8);
            app.BrowseInBtn=uibutton(app.FilesPanel,'Text','Browse','Position',[292 62 50 22],'FontSize',8,'ButtonPushedFcn',@(~,~)app.onBrowseIn());
            uilabel(app.FilesPanel,'Text','Output dir:','Position',[8 30 80 20],'FontSize',9);
            app.OutputField=uieditfield(app.FilesPanel,'text','Position',[90 30 198 22],'FontSize',8);
            app.BrowseOutBtn=uibutton(app.FilesPanel,'Text','Browse','Position',[292 30 50 22],'FontSize',8,'ButtonPushedFcn',@(~,~)app.onBrowseOut());
            y=y-114;

            % Slope
            app.SlopePanel=uipanel(app.LeftPanel,'Title','Slope & Kinematic Parameters','FontWeight','bold','FontSize',10,'Position',[6 y-132 LW-12 136]);
            sP={'Slope dip [deg]:',80,0,90,'SlopeDipSpin';'Slope dip dir [deg]:',135,0,359,'SlopeDDirSpin';'Friction angle [deg]:',35,1,89,'FrictionSpin'};
            for ii=1:3
                py=136-28-ii*34;
                uilabel(app.SlopePanel,'Text',sP{ii,1},'Position',[8 py 170 22],'FontSize',9);
                app.(sP{ii,5})=uispinner(app.SlopePanel,'Value',sP{ii,2},'Limits',[sP{ii,3} sP{ii,4}],'Step',1,'Position',[178 py 118 22],'FontSize',9);
            end
            y=y-148;

            % Algorithm
            app.AlgoPanel=uipanel(app.LeftPanel,'Title','Algorithm Parameters','FontWeight','bold','FontSize',10,'Position',[6 y-226 LW-12 230]);
            aP={'Max points [M]:',10,'MaxPtsFld';'Voxel size [m]:',0.05,'VoxelFld';'Normal k-NN:',30,'NormalKFld';'Angular eps [deg]:',8,'AngleEpsFld';'Min cluster pts:',300,'MinClFld';'Planarity min:',0.20,'PlanarMinFld';'Planar kin win [deg]:',25,'PlanarWinFld';'Toppling kin win [deg]:',30,'ToppWinFld';'MC runs:',50,'MCRunsFld'};
            for ii=1:size(aP,1)
                py=230-20-ii*23;
                uilabel(app.AlgoPanel,'Text',aP{ii,1},'Position',[8 py 165 20],'FontSize',9);
                app.(aP{ii,3})=uieditfield(app.AlgoPanel,'numeric','Value',aP{ii,2},'Position',[172 py 118 22],'FontSize',9);
            end
            y=y-242;

            % Validation sources (6 comparison sources: 2-7; source 1 = point cloud)
            %   FIX: the row pitch is now derived from the available panel height and
            %   the "Run Validation" button is given a reserved zone at the bottom, so
            %   all six rows fit. Previously rowH was hard-coded to 30 px while vsH was
            %   clamped to 208 px, which pushed the 6th row (RocFall3) down to y=8 —
            %   directly under the button (y=6..32), which is drawn afterward and hid it.
            vsH=min(y-70, 248);
            app.ValSrcPanel=uipanel(app.LeftPanel,'Title','Validation Sources (load any CSV to compare)','FontWeight','bold','FontSize',9,'Position',[6 max(y-vsH,68) LW-12 vsH]);
            srcDefs={'Manual sets (field compass):','ManualFld','ManualBtn','onBrowseManual';
                     'qFacet FM sets:','QfFMFld','QfFMBtn','onBrowseQfFM';
                     'qFacet KD sets:','QfKDFld','QfKDBtn','onBrowseQfKD';
                     'DSE sets:','DSEFld','DSEBtn','onBrowseDSE';
                     'Observed blocks (field):','ObsFld','ObsBtn','onBrowseObs';
                     'RocFall3 output:','RF3Fld','RF3Btn','onBrowseRF3'};
            nSrc   = size(srcDefs,1);
            btnH   = 26;           % Run Validation button height
            btnGap = 8;            % gap reserved between last row and the button
            topPad = 26;           % space under the panel title
            fldH   = 20;           % edit-field / browse-button height
            usable = vsH - topPad - (btnH + btnGap);      % height available for the rows
            rowH   = max(fldH+2, floor(usable / nSrc));   % adaptive pitch (never < field+2)
            topY   = vsH - topPad;
            for si=1:nSrc
                py=topY-si*rowH;
                uilabel(app.ValSrcPanel,'Text',srcDefs{si,1},'Position',[6 py+1 150 18],'FontSize',8);
                app.(srcDefs{si,2})=uieditfield(app.ValSrcPanel,'text','Position',[155 py 140 fldH],'FontSize',7);
                app.(srcDefs{si,3})=uibutton(app.ValSrcPanel,'Text','...','Position',[298 py 36 fldH],'FontSize',8,'ButtonPushedFcn',@(~,~)app.(srcDefs{si,4})());
            end
            app.RunValBtn=uibutton(app.ValSrcPanel,'Text','Run Validation (sources 2-7)','Position',[6 btnGap LW-30 btnH],'FontSize',9,'FontWeight','bold','BackgroundColor',[0.18 0.48 0.80],'FontColor','w','ButtonPushedFcn',@(~,~)app.onRunValidation());

            % RUN button
            app.RunBtn=uibutton(app.LeftPanel,'Text','  RUN GEORANK_3D  ','FontSize',14,'FontWeight','bold','BackgroundColor',[0.04 0.50 0.04],'FontColor','w','Position',[6 10 LW-12 52],'ButtonPushedFcn',@(~,~)app.onRun());
            app.StatusLbl=uilabel(app.LeftPanel,'Text','Ready','FontSize',8,'FontColor',[0.35 0.35 0.35],'Position',[6 64 LW-12 18]);

            %=============================================================
            % RIGHT PANEL — 9 result tabs
            %=============================================================
            app.RightPanel=uipanel(app.UIFigure,'Title','Results','FontWeight','bold','FontSize',11,'Position',[LW+12 4 RW H-8]);
            app.TabGroup=uitabgroup(app.RightPanel,'Position',[4 4 RW-14 H-40]);

            function tbl=makeTable(parent,pos)
                tbl=uitable(parent,'Position',pos,'RowStriping','on','FontSize',9);
                try; tbl.ColumnSortable=true; catch; end
            end

            % 1 — Set Metrics
            app.Tab_Sets=uitab(app.TabGroup,'Title',' Set Metrics ');
            app.Tbl_Sets=makeTable(app.Tab_Sets,[4 4 RW-26 H-82]);

            % 2 — Pole Plot
            app.Tab_Pole=uitab(app.TabGroup,'Title',' Pole Plot ');
            app.Ax_Pole=uiaxes(app.Tab_Pole,'Position',[30 20 RW-80 H-108]);

            % 3 — Reliability
            app.Tab_Rel=uitab(app.TabGroup,'Title',' Reliability ');
            app.Ax_Rel=uiaxes(app.Tab_Rel,'Position',[30 20 RW-80 H-108]);

            % 4 — Wedge Pairs
            app.Tab_Wedge=uitab(app.TabGroup,'Title',' Wedge Pairs ');
            app.Tbl_Wedge=makeTable(app.Tab_Wedge,[4 4 RW-26 H-82]);

            % 4b — Kinematic dispersion (population vs mean feasibility)
            app.Tab_Disp=uitab(app.TabGroup,'Title',' Kinematic Dispersion ');
            app.Tbl_DispSets=makeTable(app.Tab_Disp,[4 H/2-30 RW-26 H/2-52]);
            app.Tbl_DispWedge=makeTable(app.Tab_Disp,[4 4 RW-26 H/2-52]);

            % 5 — Multi-Source Comparison (Manual, qFacet FM, qFacet KD, DSE)
            app.Tab_MultiVal=uitab(app.TabGroup,'Title',' Set Comparison ');
            app.Tbl_MultiVal=makeTable(app.Tab_MultiVal,[4 H/2-30 RW-26 H/2-52]);
            app.Ax_DevHeat=uiaxes(app.Tab_MultiVal,'Position',[30 20 RW-80 H/2-50]);

            % 6 — Block Volume Comparison (Observed Blocks)
            app.Tab_BlockComp=uitab(app.TabGroup,'Title',' Block Comparison ');
            app.Tbl_BlockComp=makeTable(app.Tab_BlockComp,[4 H/2-30 RW-26 H/2-52]);
            app.Ax_BlockHist=uiaxes(app.Tab_BlockComp,'Position',[30 20 RW-80 H/2-50]);

            % 7 — RocFall3 Output
            app.Tab_RocFall=uitab(app.TabGroup,'Title',' RocFall3 ');
            app.Tbl_RocStats=makeTable(app.Tab_RocFall,[4 H/2-30 RW-26 H/2-52]);
            app.Ax_RocHist=uiaxes(app.Tab_RocFall,'Position',[30 20 RW-80 H/2-50]);

            % 8 — Log
            app.Tab_Log=uitab(app.TabGroup,'Title',' Run Log ');
            app.LogArea=uitextarea(app.Tab_Log,'Position',[4 4 RW-26 H-82],'FontSize',8,'Editable','off','FontName','Courier New');

            app.UIFigure.Visible='on';
        end
    end

    methods (Access = public)
        function app=GeoRank_3D_App()
            build(app); registerApp(app,app.UIFigure);
            if nargout==0; clear app; end
        end
        function delete(app); delete(app.UIFigure); end
    end
end % classdef
