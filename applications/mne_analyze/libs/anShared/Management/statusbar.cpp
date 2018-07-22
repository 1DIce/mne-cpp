//=============================================================================================================
/**
* @file     statusbar.cpp
* @author   Lars Debor <lars.debor@tu-ilmenau.de>;
*           Simon Heinke <simon.heinke@tu-ilmenau.de>;
*           Matti Hamalainen <msh@nmr.mgh.harvard.edu>
* @version  1.0
* @date     July, 2018
*
* @section  LICENSE
*
* Copyright (C) 2018, Lars Debor, Simon Heinke and Matti Hamalainen. All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that
* the following conditions are met:
*     * Redistributions of source code must retain the above copyright notice, this list of conditions and the
*       following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and
*       the following disclaimer in the documentation and/or other materials provided with the distribution.
*     * Neither the name of MNE-CPP authors nor the names of its contributors may be used
*       to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
* WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
* PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
* INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*
* @brief    StatusBar class definition.
*
*/


//*************************************************************************************************************
//=============================================================================================================
// INCLUDES
//=============================================================================================================

#include "statusbar.h"
#include <anShared/Management/communicator.h>
#include <anShared/Management/event.h>
#include <anShared/Utils/metatypes.h>


//*************************************************************************************************************
//=============================================================================================================
// QT INCLUDES
//=============================================================================================================

#include <QLabel>
#include <QDebug>


//*************************************************************************************************************
//=============================================================================================================
// USED NAMESPACES
//=============================================================================================================

using namespace ANSHAREDLIB;


//*************************************************************************************************************
//=============================================================================================================
// DEFINE GLOBAL METHODS
//=============================================================================================================


//*************************************************************************************************************
//=============================================================================================================
// DEFINE MEMBER METHODS
//=============================================================================================================

StatusBar::StatusBar(QWidget *pParent)
    : QStatusBar(pParent)
    , m_iMsgTimeout(20000)
{
    QVector<EVENT_TYPE> vSubs = {EVENT_TYPE::STATUS_BAR_MSG};
    m_pCommunicator = new Communicator(std::move(vSubs));

    this->setStyleSheet("QStatusBar { color: red}");

    connect(m_pCommunicator, &Communicator::receivedEvent,
            this, &StatusBar::onNewMessageReceived);
}


//*************************************************************************************************************

StatusBar::~StatusBar()
{
    delete m_pCommunicator;
}


//*************************************************************************************************************

void StatusBar::onNewMessageReceived(const QSharedPointer<Event> pEvent)
{
    switch(pEvent->getType()) {
        case EVENT_TYPE::STATUS_BAR_MSG: {
        if(pEvent->getData().canConvert<QString>()) {
            showMessage(pEvent->getData().toString(), m_iMsgTimeout);
            break;
        }

        }
        default:
            qDebug() << "StatusBar received an Event that is not handled by switch-cases";
            break;
    }
}


//*************************************************************************************************************