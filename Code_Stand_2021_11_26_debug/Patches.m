classdef Patches < handle
    %Pacthes defines a patch within a
    %multiple patches IGA code
    
    properties
        %geometry
        XI
        ETA
        p
        q
        KP
        w8
        
        nkn_XI
        nkn_ETA
        nekn_XI
        nekn_ETA
        XI_elem
        ETA_elem
        n
        m
        nel
        nnp
        nen
        nqp
        ndm
        ndf
        new_knots_XI
        new_knots_ETA
        IEN
        ICN
        
        
        %analysis
        
    end
    
    methods
        function obj = Patches(XI,ETA,p,q,KP,w8)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.XI = XI;
            obj.ETA = ETA;
            obj.p = p;
            obj.q = q;
            obj.KP = KP;
            obj.w8 = w8;
            
            [obj.nkn_XI,obj.nekn_XI,obj.XI_elem] = element_extractionunique(XI);
            [obj.nkn_ETA,obj.nekn_ETA,obj.ETA_elem] = element_extractionunique(ETA);
            obj.n = nkn_XI - p - 1;
            obj.m = nkn_ETA - q - 1;
            obj.nel = obj.nekn_XI * obj.nekn_ETA;
            obj.nnp = obj.n * obj.m;
            obj.nen = (p+1) * (q+1);
            obj.nqp = 4;
            obj.ndm = size(KP,2);
            obj.ndf = obj.ndm;
            obj.IEN = [];
            obj.ICN = [];
            obj.new_knots_XI = [];
            obj.new_knots_ETA = [];
            
        end
        
        
    end
end

