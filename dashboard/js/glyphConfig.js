

function glyphSettings()
{
    var out = [

        {gene:'8030451A03Rik',color: '#F0F8FF',  glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'9130024F11Rik',color: '#FAEBD7',  glyphSymbol: '.',  glyphName: 'point'},
        {gene:'9630002D21Rik',color: '#00FFFF',  glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Abca13',color: '#7FFFD4',  glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Abi3bp',color: '#F0FFFF',  glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Adam33',color: '#F5F5DC',  glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene:'Adamts18',color: '#FFE4C4',  glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'Adamts19',color: '#FFEBCD',  glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Adamtsl5',color: '#0000FF',         glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Adarb2',color: '#8A2BE2',         glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Adcyap1',color: '#A52A2A',         glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Adra1a',color: '#DEB887',         glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Adra1b',color: '#5F9EA0',         glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Aifm3',color: '#7FFF00',         glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Alkal2',color: '#D2691E',         glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Angpt1',color: '#FF7F50',         glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Anxa2',color: '#6495ED',         glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene:'Anxa4',color: '#FFF8DC',         glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Arhgap6',color: '#DC143C',         glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'Arhgef38',color: '#00FFFF',       glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Asb4',color: '#00008B',       glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Best1',color: '#008B8B',       glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Brs3',color: '#B8860B',       glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'C030005K06Rik',color: '#A9A9A9',       glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'C130060K24Rik',color: '#006400',       glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'C1ql3',color: '#BDB76B',         glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'C1ql4',color: '#8B008B',         glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'C7',color: '#556B2F',         glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'CN725425',color: '#FF8C00',         glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Cacna2d1',color: '#9932CC',         glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Calb1',color: '#8B0000',         glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Calb2',color: '#E9967A',         glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Calca',color: '#8FBC8F',         glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Calcr',color: '#483D8B',         glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Calcrl',color: '#2F4F4F',      glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Cartpt',color: '#00CED1',      glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Cbln1',color: '#9400D3',      glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Cbln2',color: '#FF1493',      glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Cbln4',color: '#00BFFF',      glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Ccdc88c',color: '#696969',        glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Cd24a',color: '#1E90FF',        glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Cd36',color: '#B22222',        glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Cdc14a',color: '#FFFAF0',        glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Cdh23',color: '#228B22',        glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Cep83',color: '#FF00FF',        glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Chat',color: '#DCDCDC',        glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Chodl',color: '#F8F8FF',        glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Chrdl1',color: '#FFD700',        glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene:'Chst9',color: '#DAA520',        glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'Cnr1',color: '#808080',        glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Cntn5',color: '#008000',         glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Col11a1',color: '#ADFF2F',         glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Col12a1',color: '#F0FFF0',         glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Col13a1',color: '#FF69B4',         glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Col14a1',color: '#CD5C5C',         glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Col15a1',color: '#4B0082',         glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Col16a1',color: '#FFFFF0',         glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Col18a1',color: '#F0E68C',         glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Col1a2',color: '#E6E6FA',         glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Col23a1',color: '#FFF0F5',         glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'Col24a1',color: '#7CFC00', glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Col25a1',color: '#FFFACD', glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Col27a1',color: '#ADD8E6', glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Col4a5',color: '#F08080', glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Col5a2',color: '#E0FFFF', glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Col6a5',color: '#FAFAD2', glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Col9a1',color: '#D3D3D3', glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene:'Cplane1',color: '#90EE90', glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Cpne4',color: '#FFB6C1', glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'Cpne7',color: '#FFA07A', glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Cpne8',color: '#20B2AA',    glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Cpne9',color: '#87CEFA',    glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Crh',color: '#778899',    glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Crhbp',color: '#B0C4DE',    glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Crhr1',color: '#FFFFE0',    glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Crim1',color: '#00FF00',    glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Cspg4',color: '#32CD32',    glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Ctss',color: '#FAF0E6',    glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Ctxn3',color: '#FF00FF',    glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Cx3cr1',color: '#800000',    glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'D130009I18Rik',color: '#66CDAA',    glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene:'Dbh',color: '#0000CD',    glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Ddc',color: '#BA55D3',    glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Dnah12',color: '#9370DB',          glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Dop1a',color: '#3CB371',          glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Ebf1',color: '#7B68EE',          glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Ebf2',color: '#00FA9A',          glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Ecel1',color: '#48D1CC',          glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Egr1',color: '#C71585',          glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Enpep',color: '#191970',          glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Epha3',color: '#F5FFFA',          glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Erbb4',color: '#FFE4E1',          glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Esr1',color: '#FFE4B5',          glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene:'Esr2',color: '#FFDEAD',          glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Etl4',color: '#000080',          glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'Eya4',color: '#FDF5E6',          glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Fbn2',color: '#808000',         glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Fbxl7',color: '#6B8E23',         glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Fgf10',color: '#FFA500',         glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Fgfr1',color: '#FF4500',         glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Foxp1',color: '#DA70D6',         glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Foxp2',color: '#EEE8AA',         glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Gabra4',color: '#98FB98',         glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Gabrg1',color: '#AFEEEE',  glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Gad1',color: '#DB7093',  glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Gad2',color: '#FFEFD5',  glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Galntl6',color: '#FFDAB9',  glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Galr1',color: '#CD853F',  glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Galr2',color: '#FFC0CB',  glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Ghr',color: '#DDA0DD',  glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Gipr',color: '#B0E0E6',  glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Glp1r',color: '#800080',  glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Gm32647',color: '#663399',  glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Gm47757',color: '#FF0000',  glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Gpc5',color: '#BC8F8F',  glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Gpr101',color: '#4169E1',  glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene:'Gpr139',color: '#8B4513',  glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'Gpr149',color: '#FA8072',  glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Gpr153',color: '#F4A460',         glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Gpr156',color: '#2E8B57',         glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Gpr176',color: '#FFF5EE',         glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Gpr39',color: '#A0522D',         glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Gpr6',color: '#C0C0C0',         glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Gpr83',color: '#87CEEB',         glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Gpr88',color: '#6A5ACD',         glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Grik1',color: '#708090',         glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Grin2c',color: '#FFFAFA',         glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene:'Grp',color: '#00FF7F',         glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Gsg1l',color: '#4682B4',         glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'Gulp1',color: '#D2B48C',       glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Hcn1',color: '#008080',       glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Hcrtr1',color: '#D8BFD8',       glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Hcrtr2',color: '#FF6347',       glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Hmcn1',color: '#40E0D0',       glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Htr2a',color: '#EE82EE',       glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'Htr2c',color: '#F5DEB3',         glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Htr7',color: '#FFFFFF',         glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Igf1',color: '#F5F5F5',         glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Ikzf1',color: '#FFFF00',         glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Il1rap',color: '#9ACD32',         glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Il1rapl2',color: '#F0F8FF',         glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Il20ra',color: '#FAEBD7',         glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Irs4',color: '#00FFFF',         glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Itga1',color: '#7FFFD4',         glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Itih5',color: '#F0FFFF',      glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Kcna1',color: '#F5F5DC',      glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Kcnh8',color: '#FFE4C4',      glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Kcnip4',color: '#FFEBCD',      glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Kcnmb2',color: '#0000FF',      glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Kcnq4',color: '#8A2BE2',        glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Klhl1',color: '#A52A2A',        glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Lama1',color: '#DEB887',        glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Lamb1',color: '#5F9EA0',        glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Lamp5',color: '#7FFF00',        glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Lepr',color: '#D2691E',        glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Lhx1',color: '#FF7F50',        glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Lhx4',color: '#6495ED',        glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Lhx9',color: '#FFF8DC',        glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene:'Lmx1a',color: '#DC143C',        glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'Lmx1b',color: '#00FFFF',        glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Maco1',color: '#00008B',         glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Maf',color: '#008B8B',         glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Mag',color: '#B8860B',         glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Maoa',color: '#A9A9A9',         glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Maob',color: '#006400',         glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Map2',color: '#BDB76B',         glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Mc4r',color: '#8B008B',         glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Mog',color: '#556B2F',         glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Mybpc1',color: '#FF8C00',         glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Myoz3',color: '#9932CC',         glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'Necab1',color: '#8B0000', glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Necab2',color: '#E9967A', glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Nefh',color: '#8FBC8F', glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Nefl',color: '#483D8B', glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Nefm',color: '#2F4F4F', glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Nfia',color: '#00CED1', glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Nfib',color: '#9400D3', glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene:'Nmu',color: '#FF1493', glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Nos1',color: '#00BFFF', glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'Npas1',color: '#696969', glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Npffr1',color: '#1E90FF',    glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Npffr2',color: '#B22222',    glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Npnt',color: '#FFFAF0',    glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Npr1',color: '#228B22',    glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Npr3',color: '#FF00FF',    glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Nps',color: '#DCDCDC',    glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Npsr1',color: '#F8F8FF',    glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Npy',color: '#FFD700',    glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Npy1r',color: '#DAA520',    glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Npy2r',color: '#808080',    glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Npy5r',color: '#008000',    glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene:'Nr2f2',color: '#ADFF2F',    glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Nr4a2',color: '#F0FFF0',    glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Nrn1',color: '#FF69B4',          glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Ntrk1',color: '#CD5C5C',          glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Nts',color: '#4B0082',          glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Onecut1',color: '#FFFFF0',          glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Onecut3',color: '#F0E68C',          glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Oprk1',color: '#E6E6FA',          glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Oprm1',color: '#FFF0F5',          glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Otof',color: '#7CFC00',          glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Otx2',color: '#FFFACD',          glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Oxtr',color: '#ADD8E6',          glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene:'P2rx7',color: '#F08080',          glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Pamr1',color: '#E0FFFF',          glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'Pax2',color: '#FAFAD2',          glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Pax3',color: '#D3D3D3',         glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Pax5',color: '#90EE90',         glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Pax6',color: '#FFB6C1',         glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Pax7',color: '#FFA07A',         glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Pde11a',color: '#20B2AA',         glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Pde7b',color: '#87CEFA',         glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Pdgfra',color: '#778899',         glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Pdyn',color: '#B0C4DE',  glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Pdzd2',color: '#FFFFE0',  glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Pdzrn3',color: '#00FF00',  glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Penk',color: '#32CD32',  glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Piezo2',color: '#FAF0E6',  glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Plk5',color: '#FF00FF',  glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Pnoc',color: '#800000',  glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Pon2',color: '#66CDAA',  glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Pon3',color: '#0000CD',  glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Postn',color: '#BA55D3',  glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Pou6f2',color: '#9370DB',  glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Prkcq',color: '#3CB371',  glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Prkd1',color: '#7B68EE',  glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene:'Prlr',color: '#00FA9A',  glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'Prox1',color: '#48D1CC',  glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Prph',color: '#C71585',         glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Prss23',color: '#191970',         glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Pth2r',color: '#F5FFFA',         glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Pvalb',color: '#FFE4E1',         glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Qrfpr',color: '#FFE4B5',         glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Rab35',color: '#FFDEAD',         glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Rai14',color: '#000080',         glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Reln',color: '#FDF5E6',         glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Rgs5',color: '#808000',         glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene:'Rln3',color: '#6B8E23',         glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Rmst',color: '#FFA500',         glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'Rnf220',color: '#FF4500',       glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Robo3',color: '#DA70D6',       glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Rorb',color: '#EEE8AA',       glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Rreb1',color: '#98FB98',       glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Rubcnl',color: '#AFEEEE',       glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Runx1',color: '#DB7093',       glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'Rxfp1',color: '#FFEFD5',         glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Rxfp2',color: '#FFDAB9',         glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Samd3',color: '#CD853F',         glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Satb1',color: '#FFC0CB',         glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Satb2',color: '#DDA0DD',         glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Scn5a',color: '#B0E0E6',         glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Scube2',color: '#800080',         glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Sema3a',color: '#663399',         glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Sema3e',color: '#FF0000',         glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Serpinb1b',color: '#BC8F8F',      glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Shox2',color: '#4169E1',      glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Six3os1',color: '#8B4513',      glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Skor2',color: '#FA8072',      glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Slc13a3',color: '#F4A460',      glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Slc17a6',color: '#2E8B57',        glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Slc17a7',color: '#FFF5EE',        glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Slc17a8',color: '#A0522D',        glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Slc18a2',color: '#C0C0C0',        glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Slc1a2',color: '#87CEEB',        glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Slc24a4',color: '#6A5ACD',        glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Slc32a1',color: '#708090',        glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Slc47a1',color: '#FFFAFA',        glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Slc4a4',color: '#00FF7F',        glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene:'Slc5a7',color: '#4682B4',        glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'Slc6a2',color: '#D2B48C',        glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Slc6a20a',color: '#008080',         glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Slc6a4',color: '#D8BFD8',         glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Slc6a5',color: '#FF6347',         glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Slco1a4',color: '#40E0D0',         glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Snap25',color: '#EE82EE',         glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Sntg2',color: '#F5DEB3',         glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Sox1ot',color: '#FFFFFF',         glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Sox2ot',color: '#F5F5F5',         glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Sox6',color: '#FFFF00',         glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Spef2',color: '#9ACD32',         glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'Sst',color: '#F0F8FF', glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'St18',color: '#FAEBD7', glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Sulf1',color: '#00FFFF', glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Synpr',color: '#7FFFD4', glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Syt10',color: '#F0FFFF', glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Tac1',color: '#F5F5DC', glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Tacr1',color: '#FFE4C4', glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene:'Tacr3',color: '#FFEBCD', glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Tcf7l2',color: '#0000FF', glyphSymbol: 'p',  glyphName: 'star5'},
        {gene:'Tfap2b',color: '#8A2BE2', glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Tfap2d',color: '#A52A2A',    glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Th',color: '#DEB887',    glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Thrb',color: '#5F9EA0',    glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Tmem178',color: '#7FFF00',    glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Tmem72',color: '#D2691E',    glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Tnc',color: '#FF7F50',    glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Tnnt1',color: '#6495ED',    glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Tnnt2',color: '#FFF8DC',    glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Tph2',color: '#DC143C',    glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Trhr',color: '#00FFFF',    glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Ttn',color: '#00008B',    glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene:'Ttr',color: '#008B8B',    glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Tut4',color: '#B8860B',    glyphSymbol: 'h',  glyphName: 'star6'},
        {gene:'Usp24',color: '#A9A9A9',          glyphSymbol: '*',  glyphName: 'asterisk'},
        {gene:'Usp43',color: '#006400',          glyphSymbol: '.',  glyphName: 'point'},
        {gene:'Vav3',color: '#BDB76B',          glyphSymbol: '+',  glyphName: 'plus'},
        {gene:'Vmn1r196',color: '#8B008B',          glyphSymbol: 'x',  glyphName: 'cross'},
        {gene:'Vmn1r206',color: '#556B2F',          glyphSymbol: 'o',  glyphName: 'circle'},
        {gene:'Vmn1r209',color: '#FF8C00',          glyphSymbol: 's',  glyphName: 'square'},
        {gene:'Vsx2',color: '#9932CC',          glyphSymbol: 'd',  glyphName: 'diamond'},
        {gene:'Vwc2',color: '#8B0000',          glyphSymbol: '>',  glyphName: 'triangleRight'},
        {gene:'Xirp2',color: '#E9967A',          glyphSymbol: 'v',  glyphName: 'triangleDown'},
        {gene:'Zeb2',color: '#8FBC8F',          glyphSymbol: '<',  glyphName: 'triangleLeft'},
        {gene:'Zfhx3',color: '#483D8B',          glyphSymbol: '^',  glyphName: 'triangleUp'},
        {gene:'Zfp804b',color: '#2F4F4F',          glyphSymbol: 'p',  glyphName: 'star5'},

        {gene: 'Generic',       taxonomy: 'generic',     glyphSymbol: 'o',  glyphName: 'circle'},

        ];

    return out
}

// //create color ramp.
// function glyphColor(y) {
//     return y === 'non_neuron' ? '#FFFFFF' : //hsv: [0 0 1]);
//         y === 'pc_or_in' ? '#407F59' :      //hsv: [.4 .5 .5]);
//             y === 'less_active' ? '#96B38F' :   //hsv: [.3 .2 .7]);
//                 y === 'pc' ? '#00FF00' :            //hsv: [1/3 1 1]);
//                     y === 'pc2' ? '#44B300' :           //hsv: [.27 1 .7]);
//                         y === 'in_general' ? '#0000FF' :    //hsv: [2/3 1 1]);
//                             y === 'sst' ? '#00B3FF' :           //hsv: [.55 1 1]);
//                                 y === 'pvalb' ? '#5C33FF' :         //hsv: [.7 .8 1]);
//                                     y === 'ngf' ? '#FF00E6' :           //hsv: [.85 1 1]);
//                                         y === 'cnr1' ? '#FF0000' :          //hsv: [ 1 1 1]);
//                                             y === 'vip' ? '#FFC700' :           //hsv: [.13 1 1]);
//                                                 y === 'cxcl14' ? '#995C00' :        //hsv: [.1 1 .6]);
//                                                     '#FD6A02';
// }