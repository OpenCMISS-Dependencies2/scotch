/* Copyright INRIA 2004
**
** This file is part of the Scotch distribution.
**
** The Scotch distribution is libre/free software; you can
** redistribute it and/or modify it under the terms of the
** GNU Lesser General Public License as published by the
** Free Software Foundation; either version 2.1 of the
** License, or (at your option) any later version.
**
** The Scotch distribution is distributed in the hope that
** it will be useful, but WITHOUT ANY WARRANTY; without even
** the implied warranty of MERCHANTABILITY or FITNESS FOR A
** PARTICULAR PURPOSE. See the GNU Lesser General Public
** License for more details.
**
** You should have received a copy of the GNU Lesser General
** Public License along with the Scotch distribution; if not,
** write to the Free Software Foundation, Inc.,
** 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
**
** $Id$
*/
/************************************************************/
/**                                                        **/
/**   NAME       : module.h                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This is the global configuration file   **/
/**                for the whole libSCOTCH library module. **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 22 jun 1998     **/
/**                                 to     13 may 1998     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to     03 oct 1998     **/
/**                # Version 3.4  : from : 01 nov 2001     **/
/**                                 to     01 nov 2001     **/
/**                # Version 4.0  : from : 12 dec 2001     **/
/**                                 to     24 nov 2005     **/
/**                                                        **/
/************************************************************/

#define MODULE_H

/*
** Debug values.
*/

#ifdef SCOTCH_DEBUG
#define SCOTCH_DEBUG_ARCH1
#define SCOTCH_DEBUG_GAIN1
#define SCOTCH_DEBUG_PARSER1
#define SCOTCH_DEBUG_BGRAPH1
#define SCOTCH_DEBUG_GEOM1
#define SCOTCH_DEBUG_GRAPH1
#define SCOTCH_DEBUG_HGRAPH1
#define SCOTCH_DEBUG_HMESH1
#define SCOTCH_DEBUG_KGRAPH1
#define SCOTCH_DEBUG_LIBRARY1
#define SCOTCH_DEBUG_MAIN1
#define SCOTCH_DEBUG_MAP1
#define SCOTCH_DEBUG_MESH1
#define SCOTCH_DEBUG_ORDER1
#define SCOTCH_DEBUG_PARSER1
#define SCOTCH_DEBUG_VGRAPH1
#define SCOTCH_DEBUG_VMESH1
#endif /* SCOTCH_DEBUG */

#ifdef SCOTCH_DEBUG_ALL
#define SCOTCH_DEBUG_ARCH1
#define SCOTCH_DEBUG_ARCH2
#define SCOTCH_DEBUG_GAIN1
#define SCOTCH_DEBUG_GAIN2
#define SCOTCH_DEBUG_PARSER1
#define SCOTCH_DEBUG_PARSER2
#define SCOTCH_DEBUG_BGRAPH1
#define SCOTCH_DEBUG_BGRAPH2
#define SCOTCH_DEBUG_GEOM1
#define SCOTCH_DEBUG_GEOM2
#define SCOTCH_DEBUG_GRAPH1
#define SCOTCH_DEBUG_GRAPH2
#define SCOTCH_DEBUG_HGRAPH1
#define SCOTCH_DEBUG_HGRAPH2
#define SCOTCH_DEBUG_HMESH1
#define SCOTCH_DEBUG_HMESH2
#define SCOTCH_DEBUG_KGRAPH1
#define SCOTCH_DEBUG_KGRAPH2
#define SCOTCH_DEBUG_LIBRARY1
#define SCOTCH_DEBUG_LIBRARY2
#define SCOTCH_DEBUG_MAIN1
#define SCOTCH_DEBUG_MAIN2
#define SCOTCH_DEBUG_MAP1
#define SCOTCH_DEBUG_MAP2
#define SCOTCH_DEBUG_MESH1
#define SCOTCH_DEBUG_MESH2
#define SCOTCH_DEBUG_ORDER1
#define SCOTCH_DEBUG_ORDER2
#define SCOTCH_DEBUG_PARSER1
#define SCOTCH_DEBUG_PARSER2
#define SCOTCH_DEBUG_VGRAPH1
#define SCOTCH_DEBUG_VGRAPH2
#define SCOTCH_DEBUG_VMESH1
#define SCOTCH_DEBUG_VMESH2
#endif /* SCOTCH_DEBUG_ALL */

/*
** Function renaming.
*/

#if ((! defined SCOTCH_COMMON_EXTERNAL) || (defined SCOTCH_COMMON_RENAME))
#define clockGet                    _SCOTCHclockGet

#define usagePrint                  _SCOTCHusagePrint

#define errorPrint                  SCOTCH_errorPrint
#define errorPrintW                 SCOTCH_errorPrintW

#define intLoad                     _SCOTCHintLoad
#define intSave                     _SCOTCHintSave
#define intAscn                     _SCOTCHintAscn
#define intPerm                     _SCOTCHintPerm
#define intRandReset                _SCOTCHintRandReset
#define intRandInit                 _SCOTCHintRandInit
/* #define intRandVal               _SCOTCHintRandVal Already a macro */
#define intSort1asc1                _SCOTCHintSort1asc1
#define intSort2asc1                _SCOTCHintSort2asc1
#define intSort2asc2                _SCOTCHintSort2asc2

#define memAllocGroup               _SCOTCHmemAllocGroup
#define memReallocGroup             _SCOTCHmemReallocGroup
#define memOffset                   _SCOTCHmemOffset
#endif /* ((! defined SCOTCH_COMMON_EXTERNAL) || (defined SCOTCH_COMMON_RENAME)) */

#ifdef SCOTCH_RENAME
#define stratdummy                  _SCOTCHstratdummy
#define stratInit                   _SCOTCHstratInit
#define stratExit                   _SCOTCHstratExit
#define stratSave                   _SCOTCHstratSave
#define stratCondEval               _SCOTCHstratCondEval
#define stratCondExit               _SCOTCHstratCondExit
#define stratCondSave               _SCOTCHstratCondSave
#define stratParserInit             _SCOTCHstratParserInit
#define stratParserInput            _SCOTCHstratParserInput
#define stratParserLex              _SCOTCHstratParserLex
#define stratParserRemain           _SCOTCHstratParserRemain
#define stratParserSelect           _SCOTCHstratParserSelect
#define stratParserError            _SCOTCHstratParserError
#define stratParserParse            _SCOTCHstratParserParse
#define stratParserParse2           _SCOTCHstratParserParse2
#define stratTestEval               _SCOTCHstratTestEval
#define stratTestExit               _SCOTCHstratTestExit
#define stratTestSave               _SCOTCHstratTestSave
#define parsermethtokentab          _SCOTCHparsermethtokentab
#define parserparamcurr             _SCOTCHparserparamcurr
#define parserstratcurr             _SCOTCHparserstratcurr
#define parserstrattab              _SCOTCHparserstrattab

#define bgraphbipartststratab       _SCOTCHbgraphbipartststratab
#define bgraphInit                  _SCOTCHbgraphInit
#define bgraphInit2                 _SCOTCHbgraphInit2
#define bgraphInit3                 _SCOTCHbgraphInit3
#define bgraphExit                  _SCOTCHbgraphExit
#define bgraphCheck                 _SCOTCHbgraphCheck
#define bgraphSwal                  _SCOTCHbgraphSwal
#define bgraphZero                  _SCOTCHbgraphZero
#define bgraphBipartEx              _SCOTCHbgraphBipartEx
#define bgraphBipartFm              _SCOTCHbgraphBipartFm
#define bgraphBipartGg              _SCOTCHbgraphBipartGg
#define bgraphBipartGp              _SCOTCHbgraphBipartGp
#define bgraphBipartMl              _SCOTCHbgraphBipartMl
#define bgraphBipartSt              _SCOTCHbgraphBipartSt
#define bgraphBipartZr              _SCOTCHbgraphBipartZr
#define bgraphStoreInit             _SCOTCHbgraphStoreInit
#define bgraphStoreExit             _SCOTCHbgraphStoreExit
#define bgraphStoreSave             _SCOTCHbgraphStoreSave
#define bgraphStoreUpdt             _SCOTCHbgraphStoreUpdt

#define gainTablInit                _SCOTCHgainTablInit
#define gainTablExit                _SCOTCHgainTablExit
#define gainTablFree                _SCOTCHgainTablFree
#define gainTablAddLin              _SCOTCHgainTablAddLin
#define gainTablAddLog              _SCOTCHgainTablAddLog
#ifdef SCOTCH_DEBUG_GAIN1                         /* If not already redefined as accelerated macro */
#define gainTablDel                 _SCOTCHgainTablDel
#endif /* SCOTCH_DEBUG_GAIN1 */
#define gainTablFrst                _SCOTCHgainTablFrst
#define gainTablNext                _SCOTCHgainTablNext
#define gainTablCheck               _SCOTCHgainTablCheck

#define geomInit                    _SCOTCHgeomInit
#define geomExit                    _SCOTCHgeomExit

#define graphInit                   _SCOTCHgraphInit
#define graphExit                   _SCOTCHgraphExit
#define graphFree                   _SCOTCHgraphFree
#define graphLoad                   _SCOTCHgraphLoad
#define graphLoad2                  _SCOTCHgraphLoad2
#define graphSave                   _SCOTCHgraphSave
#define graphBase                   _SCOTCHgraphBase
#define graphCheck                  _SCOTCHgraphCheck
#define graphCoarsen                _SCOTCHgraphCoarsen
#define graphInduceList             _SCOTCHgraphInduceList
#define graphInducePart             _SCOTCHgraphInducePart
#define graphGeomLoadChac           _SCOTCHgraphGeomLoadChac
#define graphGeomLoadHabo           _SCOTCHgraphGeomLoadHabo
#define graphGeomLoadScot           _SCOTCHgraphGeomLoadScot
#define graphGeomSaveChac           _SCOTCHgraphGeomSaveChac
#define graphGeomSaveScot           _SCOTCHgraphGeomSaveScot

#define hallOrderHdHalmd            _SCOTCHhallOrderHdHalmd
#define hallOrderHfR2hamdf4         _SCOTCHhallOrderHfR2hamdf4
#define hallOrderHxBuild            _SCOTCHhallOrderHxBuild
#define hallOrderHxTree             _SCOTCHhallOrderHxTree

#define hgraphorderststratab        _SCOTCHhgraphorderststratab
#define hgraphInit                  _SCOTCHhgraphInit
#define hgraphExit                  _SCOTCHhgraphExit
#define hgraphFree                  _SCOTCHhgraphFree
#define hgraphInduceList            _SCOTCHhgraphInduceList
#define hgraphCheck                 _SCOTCHhgraphCheck
#define hgraphOrderBl               _SCOTCHhgraphOrderBl
#define hgraphOrderCp               _SCOTCHhgraphOrderCp
#define hgraphOrderGp               _SCOTCHhgraphOrderGp
#define hgraphOrderHd               _SCOTCHhgraphOrderHd
#define hgraphOrderHf               _SCOTCHhgraphOrderHf
#define hgraphOrderHxFill           _SCOTCHhgraphOrderHxFill
#define hgraphOrderNd               _SCOTCHhgraphOrderNd
#define hgraphOrderSi               _SCOTCHhgraphOrderSi
#define hgraphOrderSt               _SCOTCHhgraphOrderSt

#define hmeshorderststratab         _SCOTCHhmeshorderststratab
#define hmeshExit                   _SCOTCHhmeshExit
#define hmeshBase                   _SCOTCHhmeshBase
#define hmeshCheck                  _SCOTCHhmeshCheck
#define hmeshInducePart             _SCOTCHhmeshInducePart
#define hmeshHgraph                 _SCOTCHhmeshHgraph
#define hmeshMesh                   _SCOTCHhmeshMesh
#define hmeshOrderBl                _SCOTCHhmeshOrderBl
#define hmeshOrderCp                _SCOTCHhmeshOrderCp
#define hmeshOrderGp                _SCOTCHhmeshOrderGp
#define hmeshOrderGr                _SCOTCHhmeshOrderGr
#define hmeshOrderHd                _SCOTCHhmeshOrderHd
#define hmeshOrderHf                _SCOTCHhmeshOrderHf
#define hmeshOrderHxFill            _SCOTCHhmeshOrderHxFill
#define hmeshOrderNd                _SCOTCHhmeshOrderNd
#define hmeshOrderSi                _SCOTCHhmeshOrderSi
#define hmeshOrderSt                _SCOTCHhmeshOrderSt

#define kgraphmapststratab          _SCOTCHkgraphmapststratab
#define kgraphInit                  _SCOTCHkgraphInit
#define kgraphExit                  _SCOTCHkgraphExit
#define kgraphMapRb                 _SCOTCHkgraphMapRb
#define kgraphMapSt                 _SCOTCHkgraphMapSt

#define mapInit                     _SCOTCHmapInit
#define mapExit                     _SCOTCHmapExit
#define mapLoad                     _SCOTCHmapLoad
#define mapSave                     _SCOTCHmapSave

#define meshInit                    _SCOTCHmeshInit
#define meshExit                    _SCOTCHmeshExit
#define meshFree                    _SCOTCHmeshFree
#define meshLoad                    _SCOTCHmeshLoad
#define meshSave                    _SCOTCHmeshSave
#define meshBase                    _SCOTCHmeshBase
#define meshGraph                   _SCOTCHmeshGraph
#define meshCoarsen                 _SCOTCHmeshCoarsen
#define meshInduceList              _SCOTCHmeshInduceList
#define meshInducePart              _SCOTCHmeshInducePart
#define meshInduceSepa              _SCOTCHmeshInduceSepa
#define meshCheck                   _SCOTCHmeshCheck
#define meshGeomLoadHabo            _SCOTCHmeshGeomLoadHabo
#define meshGeomLoadScot            _SCOTCHmeshGeomLoadScot
#define meshGeomSaveScot            _SCOTCHmeshGeomSaveScot

#define orderInit                   _SCOTCHorderInit
#define orderExit                   _SCOTCHorderExit
#define orderLoad                   _SCOTCHorderLoad
#define orderSave                   _SCOTCHorderSave
#define orderSaveMap                _SCOTCHorderSaveMap
#define orderSaveTree               _SCOTCHorderSaveTree
#define orderCheck                  _SCOTCHorderCheck
#define orderPeri                   _SCOTCHorderPeri
#define orderRang                   _SCOTCHorderRang
#define orderTree                   _SCOTCHorderTree

#define listInit                    _SCOTCHlistInit
#define listExit                    _SCOTCHlistExit
#define listAlloc                   _SCOTCHlistAlloc
#define listFree                    _SCOTCHlistFree
#define listLoad                    _SCOTCHlistLoad
#define listSave                    _SCOTCHlistSave
#define listSort                    _SCOTCHlistSort
#define listCopy                    _SCOTCHlistCopy

#define archInit                    _SCOTCHarchInit
#define archExit                    _SCOTCHarchExit
#define archFree                    _SCOTCHarchFree
#define archLoad                    _SCOTCHarchLoad
#define archSave                    _SCOTCHarchSave
/* #define archName                 _SCOTCHarchName */
#define archClass                   _SCOTCHarchClass
#define archClassTab                _SCOTCHarchClassTab
#define archDomLoad                 _SCOTCHarchDomLoad
#define archDomSave                 _SCOTCHarchDomSave
#ifdef SCOTCH_DEBUG_ARCH2                         /* If already redefined */
#define archDomNum                  _SCOTCHarchDomNum
#define archDomTerm                 _SCOTCHarchDomTerm
#define archDomSize                 _SCOTCHarchDomSize
#define archDomWght                 _SCOTCHarchDomWght
#define archDomDist                 _SCOTCHarchDomDist
#define archDomFrst                 _SCOTCHarchDomFrst
#define archDomBipart               _SCOTCHarchDomBipart
#endif /* SCOTCH_DEBUG_ARCH2 */
#define archBuild                   _SCOTCHarchBuild
#define archCmpltArchLoad           _SCOTCHarchCmpltArchLoad
#define archCmpltArchSave           _SCOTCHarchCmpltArchSave
#define archCmpltDomNum             _SCOTCHarchCmpltDomNum
#define archCmpltDomTerm            _SCOTCHarchCmpltDomTerm
#define archCmpltDomSize            _SCOTCHarchCmpltDomSize
#define archCmpltDomWght            _SCOTCHarchCmpltDomWght
#define archCmpltDomDist            _SCOTCHarchCmpltDomDist
#define archCmpltDomFrst            _SCOTCHarchCmpltDomFrst
#define archCmpltDomLoad            _SCOTCHarchCmpltDomLoad
#define archCmpltDomSave            _SCOTCHarchCmpltDomSave
#define archCmpltDomBipart          _SCOTCHarchCmpltDomBipart
#define archDecoArchBuild           _SCOTCHarchDecoArchBuild
#define archDecoArchFree            _SCOTCHarchDecoArchFree
#define archDecoArchLoad            _SCOTCHarchDecoArchLoad
#define archDecoArchSave            _SCOTCHarchDecoArchSave
#define archDecoDomNum              _SCOTCHarchDecoDomNum
#define archDecoDomTerm             _SCOTCHarchDecoDomTerm
#define archDecoDomSize             _SCOTCHarchDecoDomSize
#define archDecoDomWght             _SCOTCHarchDecoDomWght
#define archDecoDomDist             _SCOTCHarchDecoDomDist
#define archDecoDomFrst             _SCOTCHarchDecoDomFrst
#define archDecoDomLoad             _SCOTCHarchDecoDomLoad
#define archDecoDomSave             _SCOTCHarchDecoDomSave
#define archDecoDomBipart           _SCOTCHarchDecoDomBipart
#define archHcubArchLoad            _SCOTCHarchHcubArchLoad
#define archHcubArchSave            _SCOTCHarchHcubArchSave
#define archHcubDomNum              _SCOTCHarchHcubDomNum
#define archHcubDomTerm             _SCOTCHarchHcubDomTerm
#define archHcubDomSize             _SCOTCHarchHcubDomSize
#define archHcubDomWght             _SCOTCHarchHcubDomWght
#define archHcubDomDist             _SCOTCHarchHcubDomDist
#define archHcubDomFrst             _SCOTCHarchHcubDomFrst
#define archHcubDomLoad             _SCOTCHarchHcubDomLoad
#define archHcubDomSave             _SCOTCHarchHcubDomSave
#define archHcubDomBipart           _SCOTCHarchHcubDomBipart
#define archTleafArchLoad           _SCOTCHarchTleafArchLoad
#define archTleafArchSave           _SCOTCHarchTleafArchSave
#define archTleafDomNum             _SCOTCHarchTleafDomNum
#define archTleafDomTerm            _SCOTCHarchTleafDomTerm
#define archTleafDomSize            _SCOTCHarchTleafDomSize
#define archTleafDomWght            _SCOTCHarchTleafDomWght
#define archTleafDomDist            _SCOTCHarchTleafDomDist
#define archTleafDomFrst            _SCOTCHarchTleafDomFrst
#define archTleafDomLoad            _SCOTCHarchTleafDomLoad
#define archTleafDomSave            _SCOTCHarchTleafDomSave
#define archTleafDomBipart          _SCOTCHarchTleafDomBipart
#define archMesh2ArchLoad           _SCOTCHarchMesh2ArchLoad
#define archMesh2ArchSave           _SCOTCHarchMesh2ArchSave
#define archMesh2DomNum             _SCOTCHarchMesh2DomNum
#define archMesh2DomTerm            _SCOTCHarchMesh2DomTerm
#define archMesh2DomSize            _SCOTCHarchMesh2DomSize
#define archMesh2DomWght            _SCOTCHarchMesh2DomWght
#define archMesh2DomDist            _SCOTCHarchMesh2DomDist
#define archMesh2DomFrst            _SCOTCHarchMesh2DomFrst
#define archMesh2DomLoad            _SCOTCHarchMesh2DomLoad
#define archMesh2DomSave            _SCOTCHarchMesh2DomSave
#define archMesh2DomBipart          _SCOTCHarchMesh2DomBipart
#define archMesh2DomBipartO         _SCOTCHarchMesh2DomBipartO
#define archMesh2DomBipartU         _SCOTCHarchMesh2DomBipartU
#define archMesh3ArchLoad           _SCOTCHarchMesh3ArchLoad
#define archMesh3ArchSave           _SCOTCHarchMesh3ArchSave
#define archMesh3DomNum             _SCOTCHarchMesh3DomNum
#define archMesh3DomTerm            _SCOTCHarchMesh3DomTerm
#define archMesh3DomSize            _SCOTCHarchMesh3DomSize
#define archMesh3DomWght            _SCOTCHarchMesh3DomWght
#define archMesh3DomDist            _SCOTCHarchMesh3DomDist
#define archMesh3DomFrst            _SCOTCHarchMesh3DomFrst
#define archMesh3DomLoad            _SCOTCHarchMesh3DomLoad
#define archMesh3DomSave            _SCOTCHarchMesh3DomSave
#define archMesh3DomBipart          _SCOTCHarchMesh3DomBipart
#define archTermArchLoad            _SCOTCHarchTermArchLoad
#define archTermArchSave            _SCOTCHarchTermArchSave
#define archTermDomNum              _SCOTCHarchTermDomNum
#define archTermDomTerm             _SCOTCHarchTermDomTerm
#define archTermDomSize             _SCOTCHarchTermDomSize
#define archTermDomWght             _SCOTCHarchTermDomWght
#define archTermDomDist             _SCOTCHarchTermDomDist
#define archTermDomFrst             _SCOTCHarchTermDomFrst
#define archTermDomLoad             _SCOTCHarchTermDomLoad
#define archTermDomSave             _SCOTCHarchTermDomSave
#define archTermDomBipart           _SCOTCHarchTermDomBipart
#define archTorus2ArchLoad          _SCOTCHarchTorus2ArchLoad
#define archTorus2ArchSave          _SCOTCHarchTorus2ArchSave
#define archTorus2DomNum            _SCOTCHarchTorus2DomNum
#define archTorus2DomTerm           _SCOTCHarchTorus2DomTerm
#define archTorus2DomSize           _SCOTCHarchTorus2DomSize
#define archTorus2DomWght           _SCOTCHarchTorus2DomWght
#define archTorus2DomDist           _SCOTCHarchTorus2DomDist
#define archTorus2DomFrst           _SCOTCHarchTorus2DomFrst
#define archTorus2DomBipart         _SCOTCHarchTorus2DomBipart
#define archTorus2DomLoad           _SCOTCHarchTorus2DomLoad
#define archTorus2DomSave           _SCOTCHarchTorus2DomSave
#define archTorus2DomBipart         _SCOTCHarchTorus2DomBipart
#define archTorus3ArchLoad          _SCOTCHarchTorus3ArchLoad
#define archTorus3ArchSave          _SCOTCHarchTorus3ArchSave
#define archTorus3DomNum            _SCOTCHarchTorus3DomNum
#define archTorus3DomTerm           _SCOTCHarchTorus3DomTerm
#define archTorus3DomSize           _SCOTCHarchTorus3DomSize
#define archTorus3DomWght           _SCOTCHarchTorus3DomWght
#define archTorus3DomDist           _SCOTCHarchTorus3DomDist
#define archTorus3DomFrst           _SCOTCHarchTorus3DomFrst
#define archTorus3DomLoad           _SCOTCHarchTorus3DomLoad
#define archTorus3DomSave           _SCOTCHarchTorus3DomSave
#define archTorus3DomBipart         _SCOTCHarchTorus3DomBipart
#define archVcmpltArchLoad          _SCOTCHarchVcmpltArchLoad
#define archVcmpltArchSave          _SCOTCHarchVcmpltArchSave
#define archVcmpltDomNum            _SCOTCHarchVcmpltDomNum
#define archVcmpltDomTerm           _SCOTCHarchVcmpltDomTerm
#define archVcmpltDomSize           _SCOTCHarchVcmpltDomSize
#define archVcmpltDomWght           _SCOTCHarchVcmpltDomWght
#define archVcmpltDomDist           _SCOTCHarchVcmpltDomDist
#define archVcmpltDomFrst           _SCOTCHarchVcmpltDomFrst
#define archVcmpltDomBipart         _SCOTCHarchVcmpltDomBipart
#define archVcmpltDomLoad           _SCOTCHarchVcmpltDomLoad
#define archVcmpltDomSave           _SCOTCHarchVcmpltDomSave
#define archVcmpltDomBipart         _SCOTCHarchVcmpltDomBipart
#define archVhcubArchLoad           _SCOTCHarchVhcubArchLoad
#define archVhcubArchSave           _SCOTCHarchVhcubArchSave
#define archVhcubDomNum             _SCOTCHarchVhcubDomNum
#define archVhcubDomTerm            _SCOTCHarchVhcubDomTerm
#define archVhcubDomSize            _SCOTCHarchVhcubDomSize
#define archVhcubDomWght            _SCOTCHarchVhcubDomWght
#define archVhcubDomDist            _SCOTCHarchVhcubDomDist
#define archVhcubDomFrst            _SCOTCHarchVhcubDomFrst
#define archVhcubDomLoad            _SCOTCHarchVhcubDomLoad
#define archVhcubDomSave            _SCOTCHarchVhcubDomSave
#define archVhcubDomBipart          _SCOTCHarchVhcubDomBipart

#define vgraphseparateststratab     _SCOTCHvgraphseparateststratab
#define vgraphInit                  _SCOTCHvgraphInit
#define vgraphExit                  _SCOTCHvgraphExit
#define vgraphCheck                 _SCOTCHvgraphCheck
#define vgraphZero                  _SCOTCHvgraphZero
#define vgraphSeparateEs            _SCOTCHvgraphSeparateEs
#define vgraphSeparateFm            _SCOTCHvgraphSeparateFm
#define vgraphSeparateGg            _SCOTCHvgraphSeparateGg
#define vgraphSeparateGp            _SCOTCHvgraphSeparateGp
#define vgraphSeparateMl            _SCOTCHvgraphSeparateMl
#define vgraphSeparateMt            _SCOTCHvgraphSeparateMt
#define vgraphSeparateSt            _SCOTCHvgraphSeparateSt
#define vgraphSeparateTh            _SCOTCHvgraphSeparateTh
#define vgraphSeparateVw            _SCOTCHvgraphSeparateVw
#define vgraphSeparateZr            _SCOTCHvgraphSeparateZr
#define vgraphStoreInit             _SCOTCHvgraphStoreInit
#define vgraphStoreExit             _SCOTCHvgraphStoreExit
#define vgraphStoreSave             _SCOTCHvgraphStoreSave
#define vgraphStoreUpdt             _SCOTCHvgraphStoreUpdt

#define vmeshseparateststratab      _SCOTCHvmeshseparateststratab
#define vmeshExit                   _SCOTCHvmeshExit
#define vmeshCheck                  _SCOTCHvmeshCheck
#define vmeshZero                   _SCOTCHvmeshZero
#define vmeshSeparateFm             _SCOTCHvmeshSeparateFm
#define vmeshSeparateGg             _SCOTCHvmeshSeparateGg
#define vmeshSeparateGr             _SCOTCHvmeshSeparateGr
#define vmeshSeparateMl             _SCOTCHvmeshSeparateMl
#define vmeshSeparateSt             _SCOTCHvmeshSeparateSt
#define vmeshSeparateZr             _SCOTCHvmeshSeparateZr
#define vmeshStoreInit              _SCOTCHvmeshStoreInit
#define vmeshStoreExit              _SCOTCHvmeshStoreExit
#define vmeshStoreSave              _SCOTCHvmeshStoreSave
#define vmeshStoreUpdt              _SCOTCHvmeshStoreUpdt
#endif /* SCOTCH_RENAME */
