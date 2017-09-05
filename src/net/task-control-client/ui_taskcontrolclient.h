/********************************************************************************
** Form generated from reading ui file 'taskcontrolclient.ui'
**
** Created: Sat May 15 03:15:43 2010
**      by: Qt User Interface Compiler version 4.5.0
**
** WARNING! All changes made in this file will be lost when recompiling ui file!
********************************************************************************/

#ifndef UI_TASKCONTROLCLIENT_H
#define UI_TASKCONTROLCLIENT_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QComboBox>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QMainWindow>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QStatusBar>
#include <QtGui/QTabWidget>
#include <QtGui/QToolBar>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_taskControlClient
{
public:
    QWidget *centralWidget;
    QGridLayout *gridLayout;
    QTabWidget *tabWidget;
    QWidget *tab;
    QGridLayout *gridLayout_2;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label;
    QSpinBox *spinBox_RTserverIPaddress1;
    QSpinBox *spinBox_RTserverIPaddress2;
    QSpinBox *spinBox_RTserverIPaddress3;
    QSpinBox *spinBox_RTserverIPaddress4;
    QSpacerItem *horizontalSpacer_2;
    QComboBox *comboBox_serverSelection;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_2;
    QSpinBox *spinBox_RTserverPort;
    QWidget *tab_2;
    QPushButton *pushButton_checkStatus;
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *taskControlClient)
    {
        if (taskControlClient->objectName().isEmpty())
            taskControlClient->setObjectName(QString::fromUtf8("taskControlClient"));
        taskControlClient->resize(540, 208);
        centralWidget = new QWidget(taskControlClient);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        gridLayout = new QGridLayout(centralWidget);
        gridLayout->setSpacing(6);
        gridLayout->setMargin(11);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        tabWidget = new QTabWidget(centralWidget);
        tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
        tab = new QWidget();
        tab->setObjectName(QString::fromUtf8("tab"));
        gridLayout_2 = new QGridLayout(tab);
        gridLayout_2->setSpacing(6);
        gridLayout_2->setMargin(11);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        label = new QLabel(tab);
        label->setObjectName(QString::fromUtf8("label"));

        horizontalLayout_2->addWidget(label);

        spinBox_RTserverIPaddress1 = new QSpinBox(tab);
        spinBox_RTserverIPaddress1->setObjectName(QString::fromUtf8("spinBox_RTserverIPaddress1"));
        spinBox_RTserverIPaddress1->setMaximum(255);
        spinBox_RTserverIPaddress1->setValue(158);

        horizontalLayout_2->addWidget(spinBox_RTserverIPaddress1);

        spinBox_RTserverIPaddress2 = new QSpinBox(tab);
        spinBox_RTserverIPaddress2->setObjectName(QString::fromUtf8("spinBox_RTserverIPaddress2"));
        spinBox_RTserverIPaddress2->setMaximum(255);
        spinBox_RTserverIPaddress2->setValue(75);

        horizontalLayout_2->addWidget(spinBox_RTserverIPaddress2);

        spinBox_RTserverIPaddress3 = new QSpinBox(tab);
        spinBox_RTserverIPaddress3->setObjectName(QString::fromUtf8("spinBox_RTserverIPaddress3"));
        spinBox_RTserverIPaddress3->setMaximum(255);
        spinBox_RTserverIPaddress3->setValue(6);

        horizontalLayout_2->addWidget(spinBox_RTserverIPaddress3);

        spinBox_RTserverIPaddress4 = new QSpinBox(tab);
        spinBox_RTserverIPaddress4->setObjectName(QString::fromUtf8("spinBox_RTserverIPaddress4"));
        spinBox_RTserverIPaddress4->setMaximum(255);
        spinBox_RTserverIPaddress4->setValue(55);

        horizontalLayout_2->addWidget(spinBox_RTserverIPaddress4);

        horizontalSpacer_2 = new QSpacerItem(88, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_2);

        comboBox_serverSelection = new QComboBox(tab);
        comboBox_serverSelection->setObjectName(QString::fromUtf8("comboBox_serverSelection"));

        horizontalLayout_2->addWidget(comboBox_serverSelection);


        gridLayout_2->addLayout(horizontalLayout_2, 0, 0, 1, 2);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        label_2 = new QLabel(tab);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        horizontalLayout_3->addWidget(label_2);

        spinBox_RTserverPort = new QSpinBox(tab);
        spinBox_RTserverPort->setObjectName(QString::fromUtf8("spinBox_RTserverPort"));
        spinBox_RTserverPort->setMaximum(64000);
        spinBox_RTserverPort->setValue(3333);

        horizontalLayout_3->addWidget(spinBox_RTserverPort);


        gridLayout_2->addLayout(horizontalLayout_3, 1, 0, 1, 1);

        tabWidget->addTab(tab, QString());
        tab_2 = new QWidget();
        tab_2->setObjectName(QString::fromUtf8("tab_2"));
        tabWidget->addTab(tab_2, QString());

        gridLayout->addWidget(tabWidget, 0, 0, 1, 1);

        pushButton_checkStatus = new QPushButton(centralWidget);
        pushButton_checkStatus->setObjectName(QString::fromUtf8("pushButton_checkStatus"));
        pushButton_checkStatus->setEnabled(true);
        pushButton_checkStatus->setMinimumSize(QSize(0, 0));
        QPalette palette;
        QBrush brush(QColor(0, 0, 0, 255));
        brush.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::WindowText, brush);
        palette.setBrush(QPalette::Active, QPalette::Text, brush);
        QBrush brush1(QColor(0, 170, 0, 255));
        brush1.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::ButtonText, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::WindowText, brush);
        palette.setBrush(QPalette::Inactive, QPalette::Text, brush);
        palette.setBrush(QPalette::Inactive, QPalette::ButtonText, brush1);
        QBrush brush2(QColor(119, 119, 119, 255));
        brush2.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Disabled, QPalette::WindowText, brush2);
        palette.setBrush(QPalette::Disabled, QPalette::Text, brush2);
        palette.setBrush(QPalette::Disabled, QPalette::ButtonText, brush2);
        pushButton_checkStatus->setPalette(palette);
        QFont font;
        font.setPointSize(9);
        font.setBold(true);
        font.setWeight(75);
        pushButton_checkStatus->setFont(font);

        gridLayout->addWidget(pushButton_checkStatus, 1, 0, 1, 1);

        taskControlClient->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(taskControlClient);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 540, 25));
        taskControlClient->setMenuBar(menuBar);
        mainToolBar = new QToolBar(taskControlClient);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        taskControlClient->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(taskControlClient);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        taskControlClient->setStatusBar(statusBar);

        retranslateUi(taskControlClient);
        QObject::connect(comboBox_serverSelection, SIGNAL(currentIndexChanged(QString)), taskControlClient, SLOT(server_changed()));

        tabWidget->setCurrentIndex(1);


        QMetaObject::connectSlotsByName(taskControlClient);
    } // setupUi

    void retranslateUi(QMainWindow *taskControlClient)
    {
        taskControlClient->setWindowTitle(QApplication::translate("taskControlClient", "taskControlClient", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("taskControlClient", "Serwer IP address:", 0, QApplication::UnicodeUTF8));
        comboBox_serverSelection->clear();
        comboBox_serverSelection->insertItems(0, QStringList()
         << QApplication::translate("taskControlClient", "cosmos", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("taskControlClient", "szajcomp", 0, QApplication::UnicodeUTF8)
        );
        label_2->setText(QApplication::translate("taskControlClient", "Serwer port:", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab), QApplication::translate("taskControlClient", "Connection", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab_2), QApplication::translate("taskControlClient", "Processes", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        pushButton_checkStatus->setToolTip(QApplication::translate("taskControlClient", "generate trajectory", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        pushButton_checkStatus->setText(QApplication::translate("taskControlClient", "Check status", 0, QApplication::UnicodeUTF8));
        pushButton_checkStatus->setShortcut(QApplication::translate("taskControlClient", "Ctrl+C", 0, QApplication::UnicodeUTF8));
        Q_UNUSED(taskControlClient);
    } // retranslateUi

};

namespace Ui {
    class taskControlClient: public Ui_taskControlClient {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_TASKCONTROLCLIENT_H
