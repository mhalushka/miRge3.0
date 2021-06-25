from pathlib import Path


class FormatHTML:
    topHTML = """
    <!DOCTYPE HTML>
    <html>
        <head>
            <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
            <meta name="viewport" content="width=device-width, initial-scale=1">
            <title>miRge3.0</title>
            <style type="text/css">
    		button.button{
    			background: #337ab7;
    			border-radius: 12px;
    			padding: 10px 10px;
    			display: block;
    			
    			#font-weight: bold;
    			font-size: 14px;
    			color:white;
    			text-decoration: none;
    			text-shadow:0px 0px 0px #fff;
    			border:1px solid #337ab7;
    
    			box-shadow: 0px 3px 1px white inset, 0px -3px 8px white, 0px 2px 5px rgba(0, 0, 0, 0.1), 0px 8px 10px rgba(0, 0, 0, 0.1);
    			-webkit-transition:box-shadow 0.5s;
    		}
    		
    		button.button:hover{
    			box-shadow: 0px 1px 1px white inset, 0px -2px 20px white, 0px 2px 5px rgba(0, 0, 0, 0.1), 0px 8px 10px rgba(0, 0, 0, 0.1);
    		}
    		button.button:active{
    			box-shadow: 0px 1px 2px rgba(0, 0, 0, 0.5) inset, 0px -2px 20px white, 0px 1px 5px rgba(0, 0, 0, 0.1), 0px 2px 10px rgba(0, 0, 0, 0.1);
    			background:-webkit-linear-gradient(top, #d1d1d1 0%,#ECECEC 100%);
    		}
    	
    		.highcharts-figure, .highcharts-data-table table {
    			min-width: 310px; 
    			max-width: 800px;
    			margin: 1em auto;
    		}
    
    		#container {
    			height: 400px;
    		}
    
    		.highcharts-data-table table {
    			font-family: Verdana, sans-serif;
    			border-collapse: collapse;
    			border: 1px solid #EBEBEB;
    			margin: 10px auto;
    			text-align: center;
    			width: 100%;
    			max-width: 500px;
    		}
    		.highcharts-data-table caption {
    			padding: 1em 0;
    			font-size: 1.2em;
    			color: #555;
    		}
    		.highcharts-data-table th {
    			font-weight: 600;
    			padding: 0.5em;
    		}
    		.highcharts-data-table td, .highcharts-data-table th, .highcharts-data-table caption {
    			padding: 0.5em;
    		}
    		.highcharts-data-table thead tr, .highcharts-data-table tr:nth-child(even) {
    			background: #f8f8f8;
    		}
    		.highcharts-data-table tr:hover {
    			background: #f1f7ff;
       		}
            </style>
            <link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.8.2/css/all.css" integrity="sha384-oS3vJWv+0UjzBfQzYUhtDYW+Pj2yciDJxpsK1OYPAYjqT085Qq/1cq5FLXAZQ7Ay" crossorigin="anonymous">
            <link rel="stylesheet" href="https://raw.githack.com/arunhpatil/GUI/master/CSS/jquery.dataTables.min.css">
            </head>
            <body>
                <div style="background-color:#000000; height:2px; width:100%;"></div>
                <script src="https://raw.githack.com/arunhpatil/GUI/master/Code/highcharts.js"></script>
                <script src="https://raw.githack.com/arunhpatil/GUI/master/modules/exporting.js"></script>
                <script src="https://raw.githack.com/arunhpatil/GUI/master/modules/export-data.js"></script>
                <script src="https://raw.githack.com/arunhpatil/GUI/master/modules/accessibility.js"></script>
                <script src="https://raw.githack.com/arunhpatil/GUI/master/modules/histogram-bellcurve.js"></script>
                <script src="https://raw.githack.com/arunhpatil/GUI/master/modules/heatmap.js"></script>
                <script src="https://raw.githack.com/arunhpatil/GUI/master/modules/tilemap.js"></script>
                
                <p style="font-weight: bold; text-align:center; font-family: 'Source Sans Pro', sans-serif; font-size: 20px;">miRge3.0: Comprehensive analysis of small RNA sequencing Data.</p>
                
                <div style="background-color:#000000; height:2px; width:100%;"></div>
                
                <table border="0" width="100%" height="100%" style=" margin-left:auto; margin-right:auto;" >
                    <tr style=" border: solid 1px; padding: 5px; border-radius: 20px; border-color: black;">
                        <td align="center" style="padding: 5px; ">
                            <button class="button btn btn-primary btn-xs" onclick="show1();" >
                                <i class="far fa-chart-bar"></i>
                                    SmallRNA distribution
                            </button>
                        </td>
                        <td align="center" >
                            <button class="button btn btn-primary btn-xs" onclick="show2();">
                                <i class="fas fa-ruler"></i>&nbsp; Read Length &nbsp; 
                            </button>
                        </td>
                        <td align="center" >
                            <button class="button btn btn-primary btn-xs"  onclick="show3();">
                                <i class="far fa-list-alt"></i>&nbsp; isomiR results &nbsp; 
                            </button>
                        </td>
                        <td align="center" >
                            <button class="button btn btn-primary btn-xs"  onclick="show4();">
                                <i class="fas fa-chart-line"></i>&nbsp; Abundant miRNAs &nbsp; 
                            </button>
                        </td>
                        <td align="center" >
                            <button class="button btn btn-primary btn-xs"  onclick="show5();">
                                <i class="fas fa-layer-group"></i>
                                    UMI distribution
                            </button>
                        </td>
                        <td align="center" >
                            <button class="button btn btn-primary btn-xs"  onclick="show6();">
                                <i class="fas fa-envelope-open-text"></i>&nbsp;
                                    Novel miRNAs
                            </button>
                        </td>
                    </tr>
                </table>    <div style="background-color:#000000; height:2px; width:100%;"></div>
                <script>
                    function show1(){
                        document.getElementById('tab_barStack').style.display = 'block';
                        document.getElementById('tab_lengthHist').style.display = 'none';
                        document.getElementById('tab_isomir').style.display = 'none';
                        document.getElementById('tab_abundantEx').style.display = 'none';
                        document.getElementById('tab_qiagenUMI').style.display = 'none';
                        document.getElementById('tab_novelmiR').style.display = 'none';
                    }
                    function show2(){
                        document.getElementById('tab_barStack').style.display = 'none';
                        document.getElementById('tab_lengthHist').style.display = 'block';
                        document.getElementById('tab_isomir').style.display = 'none';
                        document.getElementById('tab_abundantEx').style.display = 'none';
                        document.getElementById('tab_qiagenUMI').style.display = 'none';
                        document.getElementById('tab_novelmiR').style.display = 'none';
                    }	
                    function show3(){
                        document.getElementById('tab_barStack').style.display = 'none';
                        document.getElementById('tab_lengthHist').style.display = 'none';
                        document.getElementById('tab_isomir').style.display = 'block';
                        document.getElementById('tab_abundantEx').style.display = 'none';
                        document.getElementById('tab_qiagenUMI').style.display = 'none';
                        document.getElementById('tab_novelmiR').style.display = 'none';    
                    }				
                    function show4(){
                        document.getElementById('tab_barStack').style.display = 'none';
                        document.getElementById('tab_lengthHist').style.display = 'none';
                        document.getElementById('tab_isomir').style.display = 'none';
                        document.getElementById('tab_abundantEx').style.display = 'block';
                        document.getElementById('tab_qiagenUMI').style.display = 'none';
                        document.getElementById('tab_novelmiR').style.display = 'none';
                    }							
                    function show5(){
                        document.getElementById('tab_barStack').style.display = 'none';
                        document.getElementById('tab_lengthHist').style.display = 'none';
                        document.getElementById('tab_isomir').style.display = 'none';
                        document.getElementById('tab_abundantEx').style.display = 'none';
                        document.getElementById('tab_qiagenUMI').style.display = 'block';
                        document.getElementById('tab_novelmiR').style.display = 'none';
                    }				
                    function show6(){
                        document.getElementById('tab_barStack').style.display = 'none';
                        document.getElementById('tab_lengthHist').style.display = 'none';
                        document.getElementById('tab_isomir').style.display = 'none';
                        document.getElementById('tab_abundantEx').style.display = 'none';
                        document.getElementById('tab_qiagenUMI').style.display = 'none';
                        document.getElementById('tab_novelmiR').style.display = 'block';
                    }						
                </script>
                <table border="0" style="width:100%;" >
                    <tr>
                        <td align="center" style="padding: 5px;">
                            <div id="tab_barStack">
                                <figure class="highcharts-figure" style="width:800px">
                                    <div id="smallRNADist"></div>
                                </figure>
                            </div>    """

    bottomHTML = """
                        </td>
                    </tr>
                </table>
                <script type="text/javascript" src="https://raw.githack.com/arunhpatil/GUI/master/Code/jquery-3.5.1.js"></script>
                <script type="text/javascript" src="https://raw.githack.com/arunhpatil/GUI/master/Code/jquery.dataTables.min.js"></script> 
                <script type="text/javascript" src="./index_data.js"></script>
            </body>
        </html>
    """
    
    def __init__(self, workDir):
        self.workDir = workDir

    def beginHTML(self):
        workDirFile = Path(self.workDir)/"miRge3_visualization.html"
        with open(workDirFile, "a+") as viz:
            viz.write(self.topHTML)

    def appendHTML(self, tag):
        workDirFile = Path(self.workDir)/"miRge3_visualization.html"
        with open(workDirFile, "a+") as viz:
            viz.write(tag)

    def divTabs(self, containerID):
        self.containerID = containerID
        fmt_container = """
                                    <div id='""" + self.containerID + """' ></div><br><hr><br>"""            
        self.appendHTML(fmt_container)
    
    def divTabs_splCase(self, containerID):
        self.containerID = containerID
        fmt_container = """
                                    <div id='""" + self.containerID + """' style="height: 600px; width: 700px"></div><br><hr><br>"""            
        self.appendHTML(fmt_container)

    def histReadLen(self, sampleSize):
        divHistTop = """
                            <div id="tab_lengthHist" style="display:none">
                                <figure class="highcharts-figure" style="width:800px">
                                    <p><i class="fas fa-search-plus" style="font-weight: bold; color: #337ab7;"></i>&nbsp;&nbsp;drag across x or y axis to zoom</p>"""
        divHistBottom = """
                                </figure>
                            </div>            """
        self.appendHTML(divHistTop)
        for i in range(sampleSize):
            num = i+1
            rli = "readLengthID_" + str(num)
            self.divTabs(rli)
        self.appendHTML(divHistBottom)
    
    def isomirsTab(self, sampleSize, isomir):
        workDirFile = Path(self.workDir)/"miRge3_visualization.html"
        if isomir:
            divisoTop = """
                            <div id="tab_isomir" style="display:none">
                                <figure class="highcharts-figure" style="width:800px">
                                    <div id="donutChartID"></div>"""
            self.appendHTML(divisoTop)
        
            for i in range(sampleSize+1):
                num = i+1
                if i == 0:
                    rli = "isomirDivID_" + str(num)
                    self.divTabs(rli)
                else: 
                    rli = "isomirDivID_" + str(num)
                    self.divTabs_splCase(rli)
    
            divisoBottom = """
                                </figure>
                            </div>        """
            self.appendHTML(divisoBottom)
        else: 
            isonull = """
                            <div id="tab_isomir" style="display:none">
                                <h2>No Data</h2>
                                <p>isomiR(s) was not part of the analysis</p>
                            </div>                """
            self.appendHTML(isonull)
    
    def exprTab(self, sampleSize):
        divexpTop = """
                            <div id="tab_abundantEx" style="display:none">
                                <figure class="highcharts-figure" style="width:800px">        """
        self.appendHTML(divexpTop)

        for i in range(sampleSize):
            num = i+1
            rli = "exprnDivID_" + str(num)
            self.divTabs(rli)

        divexpBottom = """
                                </figure> 
                            </div>"""
        self.appendHTML(divexpBottom)

    
    def umiTab(self, sampleSize, is_umi):
        if is_umi:
            divumiTop = """
                            <div id="tab_qiagenUMI" style="display:none">
                                <figure class="highcharts-figure" style="width:800px">        
                                    <p><i class="fas fa-search-plus" style="font-weight: bold; color: #337ab7;"></i>&nbsp;&nbsp;drag across x or y axis to zoom</p>"""
            self.appendHTML(divumiTop)
            for i in range(sampleSize):
                num = i+1
                rli = "umiDivID_" + str(num)
                self.divTabs(rli)
            divumiBottom = """
                                </figure>
                            </div>        """
            self.appendHTML(divumiBottom)
        else:
            uminull = """
                            <div id="tab_qiagenUMI" style="display:none">
                                <h2>No Data</h2>
                                <p>UMI was not part of the analysis</p>
                            </div>                """
            self.appendHTML(uminull)

    
    def novelTab(self, isnmiR):
        if isnmiR == 1: #TRUE
            nmirVar = """
                            <div id="tab_novelmiR" style="display:none; margin-right: 50px; margin-left: -100px;">
                                <figure class="highcharts-figure" style="width:800px">
                                    <table border="0" id="novelmiRData" class="display" width="80%"></table>
                                </figure>
                            </div>                """
            self.appendHTML(nmirVar)
        else: #FALSE 
            nmirVar = """
                            <div id="tab_novelmiR" style="display:none">
                                <h2>No Data</h2>
                                <p>Novel miRNAs was not part of the analysis</p>
                            </div>                """
            self.appendHTML(nmirVar)

    def closeHTML(self):
        workDirFile = Path(self.workDir)/"miRge3_visualization.html"
        with open(str(workDirFile), "a+") as viz:
            viz.write(self.bottomHTML+"\n")


class FormatJS:
    id_num=1
    hid_num = 1
    def __init__(self, workDir):
        self.workDir = workDir

    def appendJS(self, tag):
        workDirFile = Path(self.workDir)/"index_data.js"
        with open(workDirFile, "a+") as viz:
            viz.write(tag)

    def readDist(self, seriesSampleNames, seriesDataParameter):
        readVar = """
Highcharts.chart('smallRNADist', {
    chart: {
        type: 'bar'
    },
    title: {
        text: 'Read distribution'
    },
    credits: {
        enabled: false
    },
    xAxis: {
        categories: """ + seriesSampleNames + """
    },
    yAxis: {
        min: 0,
        title: {
            text: 'Total read distribution'
        }
    },
    legend: {
        reversed: true
    },
    tooltip: {
        pointFormat: '<span style="color:{series.color}">{series.name}</span>: <b>{point.y}</b> ({point.percentage:.0f}%)<br/>',
        shared: true
    },
    plotOptions: {
        series: {
            stacking: 'percent'
        }
    },
    series: """ + seriesDataParameter + """
});
        """
        self.appendJS(readVar)

    def readLenDist(self, divTagID, sampleName, hist, bins):
        readLenVar = """
var chart = Highcharts.chart('""" + divTagID + """', {
    title: {
        text: '""" + sampleName + """: Read Length Distribution'
    },
    chart: {
        marginRight: 80,
        zoomType: 'xy'
    },
    credits: {
        enabled: false
    },
    xAxis: {
        categories: """+ bins +"""
    },
    yAxis: {
        allowDecimals: false,
        title: {
            text: 'Frequency'
        }
    },
    series: [{
        name: 'Read length',
        pointWidth: 10,
        type: 'column',
        colorByPoint: false,
        data: """+ hist + """,
    }]
});
        """
        self.appendJS(readLenVar)

    def sampleUMIDist(self, divUmiTagID, sampleName, hist, bins):
        umiVar = """
var chart = Highcharts.chart('""" + divUmiTagID + """', {
    title: {
        text: '""" + sampleName + """: UMI distribution'
    },
    chart: {
        zoomType: 'xy'
    },
    credits: {
        enabled: false
    },
    xAxis: {
        categories: """ +  bins + """
    },
    yAxis: [{
        title: {
            text: ''
        }}, {
        title: {
            text: 'Frequency'
        },
    }],
    series: [ {
        type: 'column',
        pointWidth: 7,
        data: """ + hist + """,
        name: 'UMI distribution',
        yAxis: 1
    }]
});
        """
        self.appendJS(umiVar)

    def topExprnChart(self, expID, sampleName, abundantExprData):
        exprnData = """
var chart = Highcharts.chart('""" + expID + """', {
    chart: {
        type: 'tilemap',
        inverted: true,
        height: '60%',
    },

    title: {
        text: '""" + sampleName + """'
    },
	subtitle: {
        text: 'Read Per Million (RPM) values of 40 most abundant miRNAs'
    },
    credits: {
        enabled: false
    },
    
    xAxis: {
        visible: false
    },

    yAxis: {
        visible: false
    },
    exporting: {
        buttons: {
            contextButton: {
                menuItems: ["viewFullscreen", "printChart", "separator", "downloadPNG", "downloadJPEG", "downloadPDF", "downloadSVG"]
            }
        }
    },
    colorAxis: {
        dataClasses: [{
            from: 0,
            to: 2000,
            color: '#DCDCDC',
            name: '< 2k RPM' 
        }, {
            from: 2000,
            to: 5000,
            color: '#D6EAF8', //#C0C0C0
            name: '2 - 5k RPM'
        }, {
            from: 5000,
            to: 15000,
            color: '#F9EDB3',
            name: '5 - 15k RPM'
        }, {
            from: 15000,
            to: 25000,
            color: '#FFC428',
            name: '15 - 25k RPM'
        }, {
            from: 25000,
            to: 35000,
            color: '#FF7987',
            name: '25 - 35k RPM'
        },{
            from: 35000,
            color: '#FF2371',
            name: '> 35k RPM'
        }]
    },
    tooltip: {
        headerFormat: '',
        pointFormat: 'The RPM of <b> {point.name}</b> is <b>{point.value}</b>'
    },
    plotOptions: {
        series: {
            dataLabels: {
                enabled: true,
                format: '{point.hc-a2}',
                color: '#000000',
                style: {
                    textOutline: false,
                    fontSize: 9
                }
            }
        }
    },
    series: [{
        name: '',
		pointPadding: 1.2,
        data: [ """+ abundantExprData + """]
    }]
});

        """
        self.appendJS(exprnData)

    def addSerialNum(self, inString):
        id_num = self.id_num
        js_nmirVar = "[\""+str(id_num)+"\"," + inString + "],"
        FormatJS.id_num += 1
        self.appendJS(js_nmirVar)


    def openNovelmiRJSData(self):
        openmiRTag = """
var dataSet = [ """
        self.appendJS(openmiRTag)

    def closeNovelmiRJSData(self):
        closenmiRTag = """
];
 
$(document).ready(function() {
    $('#novelmiRData').DataTable( {	
        data: dataSet,
        columns: [
            { title: "id" },
            { title: "Name" },
            { title: "Probability" },
            { title: "Chr" },
            { title: "Start pos." },
            { title: "End Pos." },
            { title: "Mature miRNA sequence"},
            { title: "miRNA read Count"},
        ]
    });
});
        """
        self.appendJS(closenmiRTag)

    def donutChartJSD(self, name, val):
        dataVar = "{name:'" + name + "',y:"+ str(val)+",selected:true},"
        self.appendJS(dataVar)

    def openDoChartJSD(self):
        donutVar = """
Highcharts.chart('isomirDivID_1', {
    chart: {
        type: 'pie',
    },
    title: {
        text: 'Cumulative isomiR variant type distribution of the samples'
    },
    exporting: {
        filename: 'isomiR_variants'
    },
    credits: {
        enabled: false
    },
    plotOptions: {
        pie: {
            innerSize: 150,
            depth: 65
        }
    },
    series: [{
        name: 'Variant type',
        data: [
        """
        self.appendJS(donutVar)

    def closeDoChartJSD(self):
        closeDoVar="""
        ]
    }]
});
        """
        self.appendJS(closeDoVar)


    def openisoHmapTop(self, sampleName, iso_xAxis):
        FormatJS.hid_num += 1 
        hid_num = self.hid_num
        hidDivTagID = "isomirDivID_" + str(hid_num)
        topVar = """
Highcharts.chart('"""+ hidDivTagID +"""', {
    chart: {
        type: 'heatmap',
        marginTop: 50,
        marginBottom: 100
    },
    title: {
        text: 'Read distribution of isomiRs for the top 20 abundant miRNAs: """+ sampleName +"""'
    },
    exporting: {
        filename: '"""+ sampleName +"""_isomir'
    },
    credits: {
        enabled: false
    },
    xAxis: {
        categories: """+str(iso_xAxis)+""",
        title: "Variant type",
        labels: {
            style: {
                color: 'red',
                fontSize: 10
            }
        }
    },
    yAxis: {
        categories: ['3p: >=5', '3p: +4', '3p: +3', '3p: +2', '3p: +1', 'Canonical', '3p: -1', '3p: -2', 'iso_snv_seed','iso_snv_central_offset', 'iso_snv_central','iso_snv_central_supp', 'iso_snv','5p: -2', '5p: -1','5p: +1','5p: +2','5p: +3','5p: +4','5p: >=5'],
        title: "miRNA - isomiRs (5p -> 3p)",
    },
    colorAxis: {
        min: 0,
        stops: [
            [0, '#e3fdff'],
            [0.5, '#055dff'],
            [0.9, '#000000']
        ]
        /*
        stops: [
            [0, '#F5F5F5'],
            [0.05, '#c0f8fc'],
            [0.1, '#f7cafa'],
            [0.2, '#FAFF56'],
            [0.5, '#FC56FF'],
            [0.7, '#A798FF'],
            [0.9, '#c4463a']
        ]*/
        //minColor: '#b0ceff'
        //maxColor: '#0363ff'
        //minColor: '#FFFFFF',
        //maxColor: Highcharts.getOptions().colors[0]
    },
    legend: {
        align: 'right',
        layout: 'vertical',
        margin: 0,
        verticalAlign: 'middle',
        y: 25,
        symbolHeight: 250
    },
    tooltip: {
        formatter: function () {
            return '<b>' + this.series.xAxis.categories[this.point.x] + '</b> has log2(read counts) of <br><b>' +
            this.point.value + '</b> at <br><b>' + this.series.yAxis.categories[this.point.y] + '</b>';
        }
    },
    series: [{
        name: 'variant read counts',
        borderWidth: 1,
        data: [
        """
        self.appendJS(topVar)

    def isoHmapData(self, iso_data):
        hmapVals = """
        """+ str(iso_data) +""","""
        self.appendJS(hmapVals)

    def closeisoHmapBottom(self):
        closeVar = """
        ],
        dataLabels: {
            enabled: false,
            color: '#000000'
        }
    }],
});
        """
        self.appendJS(closeVar)

